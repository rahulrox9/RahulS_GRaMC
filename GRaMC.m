%% GRaMC – Grain Reconstruction and Multi-stage Cleaning
% =========================================================================
% Description:
%   GRaMC_extended is a grain-scale downstream pipeline for MAPClean outputs.
%   It reconstructs grains from phase-specific EBSD data, removes
%   topologically unreliable GOS = 0 grains for Anorthite, merges
%   crystallographically compatible twin-related daughter grains for
%   Anorthite, and optionally applies a final area-based population filter
%   before grain statistics.
%
%   This extended version performs:
%     1. Automatic phase splitting from MAPClean outputs
%     2. Automatic processing of all indexed phases
%     3. Initial grain reconstruction for all phases
%     4. Removal of GOS = 0 grains for Anorthite only
%     5. Twin boundary detection and twin merging for Anorthite only
%     6. Pixel-level EBSD reorientation into the merged parent frame
%     7. Grain-size histogram visualisation in log10(area) space
%     8. Area-based filtering of unwanted small-grain populations
%     9. Final phase-wise output export
%
%   Outputs:
%     - Boundary, GOS, and IPF maps at each major stage
%     - Twin-boundary classification map for Anorthite
%     - Area-distribution histogram and colour-linked grain map
%     - Stage-wise .mat checkpoints containing both EBSD and grains
%     - Organised sample/phase export folders
%
% Dependencies:
%   - MTEX Toolbox (tested on v 6.0.0)
%
% Author: Rahul Subbaraman
% Version 1.0: December 2025, 2.0: April 2026, 3.0: April 2026
% =========================================================================

clc; clear; close all;
import mtex.*;
setMTEXpref('generatingHelpMode','silent');
warning('off','all');

%% ================= USER CONFIGURATION =================

%% --- 1. Directories ---
% Sets the input, output, and checkpoint directories relative to the current
% working directory. Output folders are created automatically if missing.
dataDir       = fullfile(pwd, 'DataFiles');
exportDir     = fullfile(pwd, 'GRaMC');
checkpointDir = fullfile(pwd, 'checkpoints');

if ~exist(exportDir,'dir'), mkdir(exportDir); end
if ~exist(checkpointDir,'dir'), mkdir(checkpointDir); end

%% --- 2. Stage Control Flags ---
% Enables or disables individual pipeline stages. Disabled stages generally
% require their corresponding checkpoint files to already exist, otherwise
% the script stops to prevent silent use of incomplete data.
runPhaseSplit  = true;   % Split MAPClean output into phase-specific .ctf files
run_initial    = true;   % Initial geometric reconstruction
run_gosClean   = true;   % Remove GOS = 0 grains (Anorthite only)
run_twinMerge  = true;    % Merge twin boundaries and update EBSD pixels (Anorthite only)
run_sizeThresh = true;    % Create colour-coded grain size distribution
run_areaFilter = true;    % Area filter
run_Final      = true;    % Save final output

%% --- 3. Global Parameters ---
% Stores global reconstruction, plotting, and filtering parameters used
% across all samples and phases. Area thresholds are defined manually per
% sample-phase using keys of the form '<sampleName>_<phaseName>'.
global params
params.exportRes       = 300;          % PNG export resolution (DPI)
params.grainThreshold  = 10*degree;    % Angle threshold for grain reconstruction
params.minPixelsGrains = 5;            % Minimum pixel count to constitute a grain
params.binsPerDecade   = 10;           % Histogram bins per 10^x to 10^(x+1)

% Manual log10(area) threshold per sample-phase
params.logAreaThresh = containers.Map('KeyType','char','ValueType','double');
% ---- Sample 01a6 ----
params.logAreaThresh('01a6_Anorthite')  = 4.4;
params.logAreaThresh('01a6_Forsterite') = 4.8;
params.logAreaThresh('01a6_Diopside')   = 5.0;
% ---- Sample 4N3C ----
params.logAreaThresh('4N3C_Anorthite')  = 4.4;
params.logAreaThresh('4N3C_Forsterite') = 4.80;
% ---- Sample 024 ----
params.logAreaThresh('024_Anorthite')   = 4.4;
params.logAreaThresh('024_Forsterite')  = 4.70;
% ---- Sample 01a2 ----
params.logAreaThresh('01a2_Anorthite')  = 4.6;
params.logAreaThresh('01a2_Forsterite') = 5.0;

%% --- 4. Parallel Pool Setup (start once, reuse) ---
% Resets any stale MATLAB parallel jobs and starts a fresh parallel pool once
% at script startup. The same pool is reused during Anorthite twin-merge
% pixel reorientation to avoid repeated pool startup overhead.
cleanupParallelPool();
fprintf('\nStarting parallel pool once at script start...\n');
pool = startParallelPool();
fprintf('✔ Parallel pool active with %d workers\n', pool.NumWorkers);

%% ================= MAIN PIPELINE =================
sampleNames = {'sample1','sample2'}; % samplename must match the "filename.ctf" without the ".ctf"
fileList = cellfun(@(s) dir(fullfile(dataDir,[s '_clean.ctf'])), sampleNames, 'UniformOutput', false);
fileList = vertcat(fileList{:});

if isempty(fileList)
    cleanupParallelPool();
    error('GRaMC_extended Error: No *_clean.ctf files found in %s. Run MAPClean first.', dataDir);
end

for fi = 1:numel(fileList)
    tSample = tic;

    %% --- File Setup ---
    % Defines the active sample path, extracts the sample name from the cleaned
    % MAPClean filename, and creates the sample-level export folder.
    filePath = fullfile(dataDir, fileList(fi).name);
    [~, sampleName, ~] = fileparts(fileList(fi).name);
    sampleName = erase(sampleName, '_clean');
    sampleExportPath = fullfile(exportDir, sampleName);
    if ~exist(sampleExportPath,'dir'), mkdir(sampleExportPath); end

    fprintf('\n============================================================\n');
    fprintf('Processing Sample: %s\n', sampleName);
    fprintf('============================================================\n');

    %% --- Checkpoint Path for Phase Split ---
    % Defines the sample-level checkpoint used to record completion of the phase
    % split stage.
    splitFile = fullfile(checkpointDir, sprintf('%s_phaseSplit.mat', sampleName));

    %% --- Load Full MAPClean Output ---
    tLoad = tic;
    fprintf('Loading full cleaned EBSD data...\n');
    ebsd_full = EBSD.load(filePath, 'convertSpatial2EulerReferenceFrame').gridify;
    fprintf('✔ Full EBSD loaded: %d points (%.2f s)\n', numel(ebsd_full), toc(tLoad));

    %% --- STEP 0: Phase Split Export / Check ---
    % Detects all indexed mineral phases in the full MAPClean EBSD object and
    % exports one phase-specific CTF file per indexed phase. Pixels belonging to
    % all other phases are reset to notIndexed in each exported file. If phase
    % splitting is disabled, the code verifies that all required phase-specific
    % CTF files already exist before continuing.
    tStep0 = tic;
    fprintf('\n[STEP 0] Phase Split Export / Check...\n');

    phases = ebsd_full.mineralList;
    notIndexedId = find(strcmpi(phases,'notIndexed'));
    MinPhaseIds = setdiff(1:numel(phases), notIndexedId);
    MinPhaseNames = phases(MinPhaseIds);
    [Nrow, Ncol] = ebsd_full.size;

    if isempty(MinPhaseNames)
        warning('No indexed phases found for sample %s. Skipping sample.', sampleName);
        continue;
    end

    fprintf('Detected indexed phases for %s:\n', sampleName);
    for p = 1:numel(MinPhaseNames)
        fprintf('  -> %s\n', MinPhaseNames{p});
    end

    if runPhaseSplit
        for p = 1:numel(MinPhaseNames)
            phaseName = MinPhaseNames{p};
            phaseIdx = find(strcmp(phases, phaseName));

            ebsd_phase_export = ebsd_full;
            maskLinear = ebsd_phase_export.phaseId ~= phaseIdx;
            maskGrid = reshape(maskLinear, Nrow, Ncol);
            ebsd_phase_export(maskGrid).phaseId = notIndexedId;

            validCount = sum(ebsd_phase_export.phaseId == phaseIdx);
            fprintf('  -> %s has %d indexed points.\n', phaseName, validCount);

            outFile = fullfile(dataDir, sprintf('%s_%s.ctf', sampleName, phaseName));
            ebsd_phase_export.export(outFile);
            fprintf('  -> Saved phase: %s\n', outFile);
        end

        save(splitFile, 'filePath');
        fprintf('✔ [STEP 0] Phase split complete (%.2f s)\n', toc(tStep0));

    else
        missingFiles = {};

        for p = 1:numel(MinPhaseNames)
            phaseName = MinPhaseNames{p};
            phaseFile = fullfile(dataDir, sprintf('%s_%s.ctf', sampleName, phaseName));

            if exist(phaseFile, 'file')
                fprintf('  -> Phase file exists: %s\n', phaseFile);
            else
                missingFiles{end+1} = phaseFile; %#ok<AGROW>
            end
        end

        if ~isempty(missingFiles)
            fprintf('✖ Missing phase-split files for sample %s:\n', sampleName);
            for k = 1:numel(missingFiles)
                fprintf('    %s\n', missingFiles{k});
            end
            cleanupParallelPool();
            error('Phase files do not exist. Set runPhaseSplit = true and rerun.');
        end

        fprintf('✔ [STEP 0] All phase files exist. Continuing (%.2f s)\n', toc(tStep0));
    end

    %% --- Process each detected indexed phase automatically ---
    % Loops through every indexed phase detected in the current sample. Each
    % phase is processed independently using its own EBSD object, checkpoint
    % files, log file, and export folder. Anorthite receives additional GOS = 0
    % cleaning and twin-merging steps; all other phases follow the generic grain
    % reconstruction, visualisation, area filtering, and final export workflow.
    for p = 1:numel(MinPhaseNames)
        tPhaseTotal = tic;
        phaseName = MinPhaseNames{p};
        isAnorthite = strcmpi(phaseName, 'Anorthite');

        phaseFile = fullfile(dataDir, sprintf('%s_%s.ctf', sampleName, phaseName));
        if ~exist(phaseFile, 'file')
            fprintf('✖ Phase file missing, skipping: %s\n', phaseFile);
            continue;
        end

        phaseExportPath = fullfile(exportDir, sampleName, phaseName);
        if ~exist(phaseExportPath,'dir'), mkdir(phaseExportPath); end

        %% --- Logging ---
        % Starts a phase-specific diary file so all command-window output for the
        % current sample-phase is saved alongside the exported maps.
        diary off;
        diaryFile = fullfile(phaseExportPath, sprintf('%s_%s_logfile.txt', sampleName, phaseName));
        if exist(diaryFile, 'file'), delete(diaryFile); end
        diary(diaryFile); diary on;

        fprintf('\n------------------------------------------------------------\n');
        fprintf('Processing Sample: %s | Phase: %s\n', sampleName, phaseName);
        fprintf('------------------------------------------------------------\n');

        %% --- Phase-specific checkpoint paths ---
        % Defines all checkpoint paths for the current sample-phase. These files
        % allow individual stages to be skipped and reloaded without rerunning the
        % full pipeline.
        initialFile = fullfile(checkpointDir, sprintf('%s_%s_initialGrains.mat', sampleName, phaseName));
        gosFile     = fullfile(checkpointDir, sprintf('%s_%s_gosCleanGrains.mat', sampleName, phaseName));
        twinFile    = fullfile(checkpointDir, sprintf('%s_%s_TwinMergedGrains.mat', sampleName, phaseName));
        areaFile    = fullfile(checkpointDir, sprintf('%s_%s_areaFilteredGrains.mat', sampleName, phaseName));
        FinalFile   = fullfile(checkpointDir, sprintf('%s_%s_finalGrains.mat', sampleName, phaseName));

        %% --- Load Phase-Specific EBSD ---
        tPhaseLoad = tic;
        fprintf('Loading phase-specific EBSD: %s\n', phaseFile);
        ebsd_phase = EBSD.load(phaseFile, 'convertSpatial2EulerReferenceFrame').gridify;
        fprintf('✔ %s EBSD loaded: %d points (%.2f s)\n', phaseName, numel(ebsd_phase), toc(tPhaseLoad));

        phases_phase = ebsd_phase.mineralList;
        notIndexedId_phase = find(strcmpi(phases_phase,'notIndexed'));
        phaseId_target = find(strcmpi(phases_phase, phaseName));
        if isempty(phaseId_target)
            cleanupParallelPool();
            error('%s phase not found in %s', phaseName, phaseFile);
        end

        %% --- STEP 1: Initial Grain Reconstruction ---
        % Reconstructs grains for the current phase using the global angular
        % threshold. Grains smaller than the minimum pixel count are treated as
        % sub-resolution noise, reset to notIndexed, and grains are recalculated so
        % the EBSD and grain objects remain consistent.
        tStep1 = tic;
        if run_initial
            fprintf('\n[STEP 1] Computing initial grains for %s (Threshold: %.1f°)...\n', ...
                phaseName, params.grainThreshold/degree);

            ebsd_initial = ebsd_phase;
            [grains_initial, ebsd_initial.grainId] = calcGrains(ebsd_initial, 'angle', params.grainThreshold);
            fprintf('  -> Initial grains found: %d\n', grains_initial.length);

            smallMask = grains_initial.grainSize < params.minPixelsGrains;
            if any(smallMask)
                fprintf('  -> Removing %d small noise grains (<%d pixels).\n', ...
                    sum(smallMask), params.minPixelsGrains);

                grainIdsToRemove = grains_initial.id(smallMask);
                pixelMaskRemove = ismember(ebsd_initial.grainId, grainIdsToRemove);
                ebsd_initial(pixelMaskRemove).phaseId = notIndexedId_phase;

                [grains_initial, ebsd_initial.grainId] = calcGrains(ebsd_initial, 'angle', params.grainThreshold);
                fprintf('  -> Grain count after small-grain removal: %d\n', grains_initial.length);
            end

            save(initialFile, 'ebsd_initial', 'grains_initial', 'filePath', 'phaseName');
            fprintf('✔ Checkpoint saved: %s\n', initialFile);

            plotGrainMaps(grains_initial, phaseExportPath, sampleName, phaseName, '01_Initial');
            fprintf('✔ [STEP 1] Complete (%.2f s)\n', toc(tStep1));
        else
            checkExists(initialFile, 'STEP 1');
            S = load(initialFile, 'ebsd_initial', 'grains_initial', 'filePath', 'phaseName');
            ebsd_initial = S.ebsd_initial;
            grains_initial = S.grains_initial;
            fprintf('✔ [STEP 1] Initial checkpoint loaded: %d grains, in %.2f s\n', ...
                grains_initial.length, toc(tStep1));
        end

        %% --- STEP 2: GOS = 0 Grain Removal (Anorthite only) ---
        % Identifies Anorthite grains with zero Grain Orientation Spread, which are
        % treated in this workflow as unreliable grain objects. Pixels belonging to
        % these grains are reset to notIndexed, and grains are recalculated. For
        % non-Anorthite phases, this step is skipped and the initial reconstruction
        % is passed forward unchanged.
        tStep2 = tic;
        if run_gosClean && isAnorthite
            fprintf('\n[STEP 2] Removing GOS = 0 grains for %s...\n', phaseName);

            ebsd_gos = ebsd_initial;
            grains_gos_pre = grains_initial;

            gosMask = grains_gos_pre.GOS == 0;
            nGOS0 = sum(gosMask);
            fprintf('  -> GOS = 0 grains found: %d\n', nGOS0);

            plotGOSZeroGrains(grains_gos_pre, gosMask, phaseExportPath, sampleName, phaseName, '02a_GOS0_Grains');

            if nGOS0 > 0
                grainIdsToRemove = grains_gos_pre.id(gosMask);
                pixelMaskRemove = ismember(ebsd_gos.grainId, grainIdsToRemove);

                ebsd_gos(pixelMaskRemove).phaseId = notIndexedId_phase;
                [grains_gosclean, ebsd_gos.grainId] = calcGrains(ebsd_gos, 'angle', params.grainThreshold);
            else
                grains_gosclean = grains_gos_pre;
            end

            fprintf('  -> Grain count after GOS = 0 removal: %d\n', grains_gosclean.length);

            save(gosFile, 'ebsd_gos', 'grains_gosclean', 'gosMask', 'filePath', 'phaseName');
            fprintf('✔ Checkpoint saved: %s\n', gosFile);

            plotGrainMaps(grains_gosclean, phaseExportPath, sampleName, phaseName, '02_GOSClean');
            fprintf('✔ [STEP 2] Complete (%.2f s)\n', toc(tStep2));

        else
            if isAnorthite
                checkExists(gosFile, 'STEP 2');
                S = load(gosFile, 'ebsd_gos', 'grains_gosclean', 'gosMask', 'filePath', 'phaseName');
                ebsd_gos = S.ebsd_gos;
                grains_gosclean = S.grains_gosclean;
                gosMask = S.gosMask;
                fprintf('✔ [STEP 2] GOS-clean checkpoint loaded: %d grains, in %.2f s\n', ...
                    grains_gosclean.length, toc(tStep2));
            else
                ebsd_gos = ebsd_initial;
                grains_gosclean = grains_initial;
                gosMask = false(grains_initial.length,1);
                fprintf('✔ [STEP 2] Skipped for %s (non-Anorthite phase) (%.2f s)\n', ...
                    phaseName, toc(tStep2));
            end
        end

        %% --- STEP 3: Twin Boundary Merging + EBSD Update (Anorthite only) ---
        % Detects crystallographically compatible Anorthite twin boundaries using
        % predefined twin laws, merges linked daughter grains into parent grains,
        % and rotates daughter-pixel orientations into the selected parent reference
        % frame. The largest-area daughter is used as the reference orientation for
        % each merged parent. Non-Anorthite phases skip this stage unchanged.
        tStep3 = tic;
        if run_twinMerge && isAnorthite
            fprintf('\n[STEP 3] Computing Twin Merges for %s...\n', phaseName);

            ebsd_twin = ebsd_gos;
            grains_preTwin = grains_gosclean;
            nGrainsBefore = grains_preTwin.length;
            fprintf('  -> Grain count before twin merge: %d\n', nGrainsBefore);

            gB = grains_preTwin.boundary;
            cs = ebsd_twin(phaseName).CS;

            % --- Twin-law definitions ---
            twinLaws = {
                {'Albite',        orientation.byAxisAngle(vector3d(Miller(0,1,0,cs)), 180*degree, cs), 5*degree}
                {'Pericline',     orientation.byAxisAngle(vector3d(Miller(0,1,0,cs,'uvw')), 180*degree, cs), 5*degree}
                {'Carlsbad',      orientation.byAxisAngle(vector3d(Miller(0,0,1,cs,'uvw')), 180*degree, cs), 5*degree}
                {'Manebach',      orientation(reflection(Miller(0,0,1,cs))), 5*degree}
                {'Baveno-Right',  orientation(reflection(Miller(0,2,1,cs))), 5*degree}
                {'Baveno-Left',   orientation(reflection(Miller(0,-2,1,cs))), 5*degree}
                {'Ala_A',         orientation.byAxisAngle(vector3d(Miller(1,0,0,cs,'uvw')), 180*degree, cs), 5*degree}
                {'Ala_B',         orientation(reflection(Miller(1,0,0,cs))), 5*degree}
            };
            
            % --- Twin boundary detection ---
            tTwinDetect = tic;
            twinId = zeros(size(gB));
            for i = 1:length(twinLaws)
                law = twinLaws{i};
                isTwin = angle(gB.misorientation, law{2}) < law{3};
                newTwins = isTwin & (twinId == 0);
                twinId(newTwins) = i;
                fprintf('  -> Detected %-14s boundaries: %d\n', law{1}, sum(newTwins));
            end
            fprintf('  -> Twin detection complete (%.2f s)\n', toc(tTwinDetect));
            
            % --- Grain-object merge ---
            plotTwinBoundaryMap(grains_preTwin, gB, twinId, twinLaws, phaseExportPath, sampleName, phaseName, '03a_TwinClassification');
            twinBoundaries = gB(twinId > 0);
            gidpair = unique(sort(twinBoundaries.grainId, 2), 'rows');
            if isempty(gidpair)
                fprintf('  -> No twin-linked grain pairs found. Skipping merge.\n');
                grains_twin = grains_preTwin;
                fprintf('  -> Grain count after twin merge: %d\n', grains_twin.length);
                save(twinFile, 'ebsd_twin', 'grains_twin', 'twinId', 'filePath', 'phaseName');
                fprintf('✔ Checkpoint saved: %s\n', twinFile);
                plotGrainMaps(grains_twin, phaseExportPath, sampleName, phaseName, '03_TwinMerge');
                fprintf('✔ [STEP 3] Complete (%.2f s)\n', toc(tStep3));
            else
                tMerge = tic;
                [grains_twin_raw, parentId] = mergeGrainsAndComputeGOS(grains_preTwin, gidpair);
                fprintf('  -> Grain-object merge complete (%.2f s)\n', toc(tMerge));

                oldGrainId = ebsd_twin.grainId;
                validPix = oldGrainId > 0 & oldGrainId <= numel(parentId);
                validPixIdx = find(validPix);
                pixByPreId = accumarray(oldGrainId(validPix), validPixIdx, [numel(parentId), 1], @(x){x});

                postIdPerPreId = parentId(:);
                uniquePostIds = unique(postIdPerPreId);
                uniquePostIds = uniquePostIds(uniquePostIds > 0);
                nPost = numel(uniquePostIds);

                refDaughterByPost = zeros(max(uniquePostIds),1);
                needsWorkPost = false(nPost,1);

                [~, sortIdx] = sort(grains_preTwin.area, 'descend');
                sortedPostId = postIdPerPreId(sortIdx);
                [uniquePostSorted, firstIdx] = unique(sortedPostId, 'stable');
                refDaughterList = sortIdx(firstIdx);
                refDaughterByPost(uniquePostSorted) = refDaughterList;

                for pp = 1:nPost
                    thisPostId = uniquePostIds(pp);
                    daughters = find(postIdPerPreId == thisPostId);
                    if numel(daughters) > 1
                        needsWorkPost(pp) = true;
                    end
                end

                workPostIds = uniquePostIds(needsWorkPost);
                nWorkPost = numel(workPostIds);

                fprintf('  -> Merged parent grains: %d\n', grains_twin_raw.length);
                fprintf('  -> Parent grains with rotating daughters: %d\n', nWorkPost);
                
                % --- Parallel pixel reorientation ---
                pool = gcp('nocreate');
                if isempty(pool)
                    cleanupParallelPool();
                    error('Parallel pool is not active. Pool should have been started at script start.');
                end
                fprintf('  -> Reusing parallel pool with %d workers\n', pool.NumWorkers);

                rotatedPixelTotal = 0;
                preMeanOri = grains_preTwin.meanOrientation;

                if nWorkPost > 0
                    tRot = tic;
                    nChunks = min(20, max(1, ceil(nWorkPost / max(pool.NumWorkers,1))));
                    chunkEdges = round(linspace(0, nWorkPost, nChunks + 1));

                    for cc = 1:nChunks
                        i1 = chunkEdges(cc) + 1;
                        i2 = chunkEdges(cc + 1);

                        if i1 > i2
                            continue;
                        end

                        thisPostChunk = workPostIds(i1:i2);
                        nThis = numel(thisPostChunk);

                        chunkPixIdx = cell(nThis,1);
                        chunkQuat   = cell(nThis,1);
                        chunkPostId = cell(nThis,1);

                        refDaughterVec   = zeros(nThis,1);
                        rotDaughtersCell = cell(nThis,1);
                        qRefCell         = cell(nThis,1);
                        qPixCell         = cell(nThis,1);
                        pixIdxCell       = cell(nThis,1);

                        % Build compact per-iteration inputs on client
                        for kk = 1:nThis
                            thisPostId = thisPostChunk(kk);
                            daughters = find(postIdPerPreId == thisPostId);
                            refDaughter = refDaughterByPost(thisPostId);
                            rotDaughters = daughters(daughters ~= refDaughter);

                            refDaughterVec(kk) = refDaughter;
                            rotDaughtersCell{kk} = rotDaughters;
                            qRefCell{kk} = quaternion(preMeanOri(refDaughter));

                            thisPixIdxCell = cell(numel(rotDaughters),1);
                            thisQPixCell   = cell(numel(rotDaughters),1);

                            for dd = 1:numel(rotDaughters)
                                gid = rotDaughters(dd);
                                pixIdx = pixByPreId{gid};
                                thisPixIdxCell{dd} = pixIdx;

                                if ~isempty(pixIdx)
                                    thisQPixCell{dd} = quaternion(ebsd_twin(pixIdx).orientations);
                                else
                                    thisQPixCell{dd} = quaternion;
                                end
                            end

                            pixIdxCell{kk} = thisPixIdxCell;
                            qPixCell{kk}   = thisQPixCell;
                        end

                        parfor kk = 1:nThis
                            thisPostId = thisPostChunk(kk);
                            refDaughter = refDaughterVec(kk);
                            rotDaughters = rotDaughtersCell{kk};
                            qRef = qRefCell{kk};
                            thisPixIdxCell = pixIdxCell{kk};
                            thisQPixCell   = qPixCell{kk};

                            pixIdxAll = [];
                            quatAll = zeros(0,4);
                            postIdAll = [];

                            for dd = 1:numel(rotDaughters)
                                gid = rotDaughters(dd);
                                pixIdx = thisPixIdxCell{dd};

                                if isempty(pixIdx)
                                    continue;
                                end

                                qD = quaternion(preMeanOri(gid));
                                rotOp = qRef * inv(qD);
                                qPix = thisQPixCell{dd};
                                qRot = rotOp * qPix;

                                pixIdxAll = [pixIdxAll; pixIdx(:)]; %#ok<AGROW>
                                quatAll = [quatAll; [qRot.a(:), qRot.b(:), qRot.c(:), qRot.d(:)]]; %#ok<AGROW>
                                postIdAll = [postIdAll; repmat(thisPostId, numel(pixIdx), 1)]; %#ok<AGROW>
                            end

                            pixIdxRef = pixByPreId{refDaughter};
                            if ~isempty(pixIdxRef)
                                pixIdxAll = [pixIdxAll; pixIdxRef(:)]; %#ok<AGROW>
                                quatAll = [quatAll; zeros(numel(pixIdxRef),4)]; %#ok<AGROW>
                                postIdAll = [postIdAll; repmat(thisPostId, numel(pixIdxRef), 1)]; %#ok<AGROW>
                            end

                            chunkPixIdx{kk} = pixIdxAll;
                            chunkQuat{kk}   = quatAll;
                            chunkPostId{kk} = postIdAll;
                        end

                        for kk = 1:nThis
                            pixIdxAll = chunkPixIdx{kk};
                            quatAll   = chunkQuat{kk};
                            postIdAll = chunkPostId{kk};

                            if isempty(pixIdxAll)
                                continue;
                            end

                            rotMask = any(quatAll ~= 0, 2);
                            if any(rotMask)
                                qObj = quaternion(quatAll(rotMask,1), quatAll(rotMask,2), quatAll(rotMask,3), quatAll(rotMask,4));
                                ebsd_twin(pixIdxAll(rotMask)).orientations = orientation(qObj, cs);
                            end

                            ebsd_twin.grainId(pixIdxAll) = postIdAll;
                            rotatedPixelTotal = rotatedPixelTotal + sum(rotMask);
                        end

                        pct = floor(100 * i2 / max(nWorkPost,1));
                        fprintf('      %3d%% | parents = %d/%d | rotated pixels = %d | elapsed = %.1f s\n', ...
                            pct, i2, nWorkPost, rotatedPixelTotal, toc(tRot));
                    end

                    fprintf('  -> Pixel rotation complete (%.2f s)\n', toc(tRot));
                end
                
                % --- Recalculate grains after EBSD update ---
                tRecalc = tic;
                [grains_twin, ebsd_twin.grainId] = calcGrains(ebsd_twin, 'angle', params.grainThreshold);
                fprintf('  -> Recalculated grains from updated EBSD (%.2f s)\n', toc(tRecalc));
                fprintf('  -> Grain count after twin merge: %d\n', grains_twin.length);

                save(twinFile, 'ebsd_twin', 'grains_twin', 'twinId', 'filePath', 'phaseName');
                fprintf('✔ Checkpoint saved: %s\n', twinFile);

                plotGrainMaps(grains_twin, phaseExportPath, sampleName, phaseName, '03_TwinMerge');
                fprintf('✔ [STEP 3] Complete (%.2f s)\n', toc(tStep3));
            end

        else
            if isAnorthite
                checkExists(twinFile, 'STEP 3');
                S = load(twinFile, 'ebsd_twin', 'grains_twin', 'twinId', 'filePath', 'phaseName');
                ebsd_twin = S.ebsd_twin;
                grains_twin = S.grains_twin;
                twinId = S.twinId;
                fprintf('✔ [STEP 3] Twin merge checkpoint loaded: %d grains, in %.2f s\n', ...
                    grains_twin.length, toc(tStep3));
            else
                ebsd_twin = ebsd_gos;
                grains_twin = grains_gosclean;
                twinId = [];
                fprintf('✔ [STEP 3] Skipped for %s (non-Anorthite phase) (%.2f s)\n', ...
                    phaseName, toc(tStep3));
            end
        end

        %% --- STEP 4: log10(Area) Histogram + Matching Grain-Bin Map ---
        % Builds a diagnostic grain-size distribution in log10(area) space and
        % exports a matching spatial grain map using the same colour bins. This helps
        % connect histogram populations to their physical locations in the EBSD map
        % before applying the manual area threshold.
                tStep4 = tic;
        if run_sizeThresh
            fprintf('\n[STEP 4] log10(area) histogram + matching grain-bin map for %s...\n', phaseName);

            areas = grains_twin.area;
            logAreas = log10(areas);

            binsPerDecade = params.binsPerDecade;
            logMin = floor(min(logAreas));
            logMax = ceil(max(logAreas));
            edges = logMin : (1 / binsPerDecade) : logMax;

            if edges(end) < max(logAreas)
                edges(end+1) = edges(end) + 1 / binsPerDecade;
            end

            nbins = numel(edges) - 1;
            binId = discretize(logAreas, edges);
            cmap = plasma(nbins);

            fHist = figure('Visible','off');
            hold on;
            counts = histcounts(logAreas, edges);
            for b = 1:nbins
                xc = (edges(b) + edges(b+1)) / 2;
                bw = (edges(b+1) - edges(b));
                bar(xc, counts(b), bw, ...
                    'FaceColor', cmap(b,:), ...
                    'EdgeColor', 'k');
            end
            xlabel('log10(grain area)');
            ylabel('Count');
            grid on;
            hold off;
            savePNG(fHist, sprintf('%s_%s_04_LogAreaHistogram', sampleName, phaseName), phaseExportPath);

            grainColors = zeros(grains_twin.length, 3);
            for b = 1:nbins
                thisMask = (binId == b);
                if any(thisMask)
                    grainColors(thisMask, :) = repmat(cmap(b,:), sum(thisMask), 1);
                end
            end

            fMap = figure('Visible','off');
            plot(grains_twin, grainColors);
            hold on;
            plot(grains_twin.boundary, 'lineColor', 'k', 'lineWidth', 0.5);
            axis equal tight;
            hold off;

            colormap(cmap);
            cb = colorbar;
            cb.Ticks = linspace(1/(2*nbins), 1 - 1/(2*nbins), nbins);
            cb.TickLabels = arrayfun(@(i) sprintf('%.2f–%.2f', edges(i), edges(i+1)), 1:nbins, 'UniformOutput', false);
            cb.Label.String = 'log10(grain area) bins';
            cb.FontSize = 6;
            cb.Label.FontSize = 6;

            savePNG(fMap, sprintf('%s_%s_04_LogAreaBinMap', sampleName, phaseName), phaseExportPath);

            fprintf('  -> Grain area stats:\n');
            fprintf('     min area        = %.3f\n', min(areas));
            fprintf('     max area        = %.3f\n', max(areas));
            fprintf('     median area     = %.3f\n', median(areas));
            fprintf('     min log10(area) = %.3f\n', min(logAreas));
            fprintf('     max log10(area) = %.3f\n', max(logAreas));
            fprintf('     bins/decade     = %d\n', binsPerDecade);

            fprintf('✔ [STEP 4] Complete (%.2f s)\n', toc(tStep4));
        else
            fprintf('✔ [STEP 4] Skipped (%.2f s)\n', toc(tStep4));
        end

        %% --- STEP 5: Area-Based Grain Filtering ---
        % Applies the manually defined log10(area) threshold for the current
        % sample-phase. Grains below the threshold are removed by resetting their
        % pixels to notIndexed, then grains are recalculated so EBSD and grain
        % objects remain consistent. A threshold of zero skips filtering.
        tStep5 = tic;
        if run_areaFilter
            threshKey = sprintf('%s_%s', sampleName, phaseName);

            if ~isKey(params.logAreaThresh, threshKey)
                cleanupParallelPool();
                error('No log10(area) threshold defined for sample-phase %s', threshKey);
            end

            logAreaThresh = params.logAreaThresh(threshKey);

            if logAreaThresh == 0
                fprintf('\n[STEP 5] Area filter skipped because log10(area) threshold = 0 for %s\n', threshKey);
                ebsd_area = ebsd_twin;
                grains_area = grains_twin;
                fprintf('✔ [STEP 5] Skipped (%.2f s)\n', toc(tStep5));
            else
                fprintf('\n[STEP 5] Area-based grain filtering for %s...\n', phaseName);

                areaThresh = 10^logAreaThresh;

                fprintf('  -> log10(area) threshold for %s = %.6f\n', threshKey, logAreaThresh);
                fprintf('  -> area threshold for %s = %.6f\n', threshKey, areaThresh);

                ebsd_area = ebsd_twin;
                grains_area_pre = grains_twin;

                fprintf('  -> Grain count before area filter: %d\n', grains_area_pre.length);

                logAreas = log10(grains_area_pre.area);
                areaMask = logAreas < logAreaThresh;

                nRemove = sum(areaMask);
                fprintf('  -> Grains below threshold: %d\n', nRemove);

                if nRemove > 0
                    grainIdsToRemove = grains_area_pre.id(areaMask);
                    pixelMaskRemove = ismember(ebsd_area.grainId, grainIdsToRemove);

                    phases_area = ebsd_area.mineralList;
                    notIndexedId_area = find(strcmpi(phases_area,'notIndexed'));

                    ebsd_area(pixelMaskRemove).phaseId = notIndexedId_area;
                    [grains_area, ebsd_area.grainId] = calcGrains(ebsd_area, 'angle', params.grainThreshold);
                else
                    grains_area = grains_area_pre;
                end

                fprintf('  -> Grain count after area filter: %d\n', grains_area.length);

                save(areaFile, 'ebsd_area', 'grains_area', 'logAreaThresh', 'areaThresh', 'filePath', 'phaseName');
                fprintf('✔ Checkpoint saved: %s\n', areaFile);

                plotGrainMaps(grains_area, phaseExportPath, sampleName, phaseName, '05_AreaFiltered');

                fprintf('✔ [STEP 5] Complete (%.2f s)\n', toc(tStep5));
            end
        else
            checkExists(areaFile, 'STEP 5');
            S = load(areaFile, 'ebsd_area', 'grains_area', 'logAreaThresh', 'areaThresh', 'filePath', 'phaseName');
            ebsd_area = S.ebsd_area;
            grains_area = S.grains_area;
            logAreaThresh = S.logAreaThresh;
            areaThresh = S.areaThresh;
            fprintf('✔ [STEP 5] Loaded checkpoint: %d grains.\n', grains_area.length);
            fprintf('✔ [STEP 5] Complete (%.2f s)\n', toc(tStep5));
        end

        %% --- STEP 6: Final Output ---
        % Selects the most advanced available EBSD/grain state for the current
        % sample-phase, saves it as the final checkpoint, and exports a cleaned
        % phase-specific CTF file for downstream analysis.
        tStep6 = tic;
        if run_Final
            if exist('ebsd_area', 'var') && exist('grains_area', 'var')
                ebsd_final = ebsd_area;
                finalGrains = grains_area;
            elseif exist('ebsd_twin', 'var') && exist('grains_twin', 'var')
                ebsd_final = ebsd_twin;
                finalGrains = grains_twin;
            elseif exist('ebsd_gos', 'var') && exist('grains_gosclean', 'var')
                ebsd_final = ebsd_gos;
                finalGrains = grains_gosclean;
            else
                ebsd_final = ebsd_initial;
                finalGrains = grains_initial;
            end

            save(FinalFile, 'ebsd_final', 'finalGrains', 'filePath', 'phaseName');
            fprintf('\n[FINAL] Final grains saved to: %s\n', FinalFile);

            outCTF = fullfile(dataDir, sprintf('%s_%s_clean.ctf', sampleName, phaseName));
            ebsd_final.export(outCTF);
            fprintf('✔ Exported cleaned CTF: %s\n', outCTF);

            fprintf('✔ [STEP 6] Complete (%.2f s)\n', toc(tStep6));
        else
            fprintf('✔ [STEP 6] Skipped (%.2f s)\n', toc(tStep6));
        end

        fprintf('✔ Phase complete: %s | %s (%.2f s)\n', sampleName, phaseName, toc(tPhaseTotal));
        diary off;
    end

    fprintf('\n✔ Sample complete: %s (%.2f s)\n', sampleName, toc(tSample));
end

fprintf('\nCleaning up parallel pool...\n');
cleanupParallelPool();
fprintf('✔ Parallel cleanup complete.\n');

%% ================= HELPER FUNCTIONS =================
% Local utility functions used by the main GRaMC_extended pipeline for
% checkpoint validation, grain merging, plotting, colormap generation, PNG
% export, and parallel-pool management.

function checkExists(fileP, stepName)
% CHECKEXISTS Verifies that a required checkpoint file exists.
%
% Inputs:
%   fileP    - Full path to the checkpoint file.
%   stepName - Name of the pipeline step requesting the checkpoint.
%
% Throws:
%   Error if the checkpoint file is missing.
    if ~exist(fileP, 'file')
        error('%s disabled but checkpoint file missing: %s', stepName, fileP);
    end
end

function [mergedGrains, parentId] = mergeGrainsAndComputeGOS(grains, gidpair)
% MERGEGRAINSANDCOMPUTEGOS Merges twin-linked grains and updates GOS.
%
% Inputs:
%   grains  - Input MTEX grain object before merging.
%   gidpair - Nx2 array of grain-ID pairs connected by twin boundaries.
%
% Outputs:
%   mergedGrains - MTEX grain object after merging linked daughter grains.
%   parentId     - Mapping from original grain IDs to merged parent IDs.
%
% Notes:
%   The largest-area daughter grain supplies the parent mean orientation.
%   Parent GOS is recalculated as an area-weighted daughter-grain average.
    [mergedGrains, parentId] = merge(grains, gidpair);
    [~, sortIdx] = sort(grains.area, 'descend');
    sortedParentId = parentId(sortIdx);
    [uParents, uIdx] = unique(sortedParentId, 'stable');
    bestChildIdx = sortIdx(uIdx);
    newOri = mergedGrains.meanOrientation;
    newOri(uParents) = grains.meanOrientation(bestChildIdx);
    mergedGrains.meanOrientation = newOri;
    weightedSum = accumarray(parentId, grains.GOS .* grains.area);
    totalArea   = accumarray(parentId, grains.area);
    mergedGrains.prop.GOS = weightedSum ./ totalArea;
end

function plotGrainMaps(grains, exportPath, sampleName, phaseName, keyWord)
% PLOTGRAINMAPS Exports standard grain maps for one processing stage.
%
% Inputs:
%   grains     - MTEX grain object to plot.
%   exportPath - Folder where PNG files are saved.
%   sampleName - Current sample name.
%   phaseName  - Current mineral phase name.
%   keyWord    - Stage label used in output filenames.
%
% Outputs:
%   Saves boundary, GOS, and IPF PNG maps to exportPath.
    % Grain Boundary
    fBound = figure('Visible','off');
    plot(grains.boundary, 'lineColor','k');
    axis equal tight;
    savePNG(fBound, sprintf('%s_%s_%s_Boundaries', sampleName, phaseName, keyWord), exportPath);

    % GOS
    fGOS = figure('Visible','off');
    plot(grains, grains.prop.GOS./degree, 'backgroundColor','k');
    mtexColorbar;
    colormap(gca, flipud(magma));
    setColorRange([0 max(grains.prop.GOS./degree)]);
    axis equal tight;
    savePNG(fGOS, sprintf('%s_%s_%s_GOS', sampleName, phaseName, keyWord), exportPath);

    % IPF
    fIPF = figure('Visible','off');
    ipfKey = ipfColorKey(grains);
    ipfKey.inversePoleFigureDirection = vector3d.Z;
    colors = ipfKey.orientation2color(grains.meanOrientation);
    plot(grains, colors);
    hold on;
    plot(grains.boundary, 'lineColor','k', 'lineWidth',1);
    hold off;
    axis equal tight;
    savePNG(fIPF, sprintf('%s_%s_%s_IPF', sampleName, phaseName, keyWord), exportPath);
end

function plotGOSZeroGrains(grains, gosMask, exportPath, sampleName, phaseName, keyWord)
% PLOTGOSZEROGRAINS Exports a diagnostic map of GOS = 0 grains.
%
% Inputs:
%   grains     - MTEX grain object before GOS cleaning.
%   gosMask    - Logical mask identifying grains with GOS = 0.
%   exportPath - Folder where the PNG file is saved.
%   sampleName - Current sample name.
%   phaseName  - Current mineral phase name.
%   keyWord    - Stage label used in output filename.
%
% Outputs:
%   Saves a PNG map highlighting GOS = 0 grains in red.
    fGOS0 = figure('Visible','off');
    hold on;
    if any(gosMask)
        plot(grains(gosMask), 'FaceColor', [1 0 0], 'EdgeColor', 'none');
    end
    plot(grains.boundary, 'lineColor', 'k', 'lineWidth', 0.5);
    axis equal tight;
    hold off;
    savePNG(fGOS0, sprintf('%s_%s_%s', sampleName, phaseName, keyWord), exportPath);
end

function plotTwinBoundaryMap(~, gB, twinId, twinLaws, exportPath, sampleName, phaseName, keyWord)
% PLOTTWINBOUNDARYMAP Exports a classified twin-boundary map.
%
% Inputs:
%   gB         - MTEX grain-boundary object.
%   twinId     - Integer boundary classification array; zero means non-twin.
%   twinLaws   - Cell array containing twin-law names and operators.
%   exportPath - Folder where the PNG file is saved.
%   sampleName - Current sample name.
%   phaseName  - Current mineral phase name.
%   keyWord    - Stage label used in output filename.
%
% Outputs:
%   Saves a PNG map showing non-twin boundaries and each detected twin-law
%   boundary class using separate colours.
    fTwins = figure('Visible','off');
    plot(gB, 'lineColor', 'k');
    hold on;
    colors = lines(length(twinLaws));
    for i = 1:length(twinLaws)
        thisTwin = gB(twinId == i);
        if isempty(thisTwin)
            continue;
        end
        plot(thisTwin, 'lineColor', colors(i,:), 'lineWidth', 2);
    end
    h = [];
    labels = {};
    h(end+1) = plot(nan, nan, '-', 'Color', 'k', 'LineWidth', 0.5);
    labels{end+1} = 'Other grain boundaries';
    for i = 1:length(twinLaws)
        h(end+1) = plot(nan, nan, '-', 'Color', colors(i,:), 'LineWidth', 2);
        labels{end+1} = twinLaws{i}{1};
    end
    legend(h, labels, 'Location', 'eastoutside');
    axis equal tight;
    savePNG(fTwins, sprintf('%s_%s_%s_TwinBoundaries', sampleName, phaseName, keyWord), exportPath);
end

function savePNG(figHandle, filenameStem, exportPath)
% SAVEPNG Saves and closes a figure using the global export resolution.
%
% Inputs:
%   figHandle    - MATLAB figure handle to export.
%   filenameStem - Output filename without extension.
%   exportPath   - Folder where the PNG file is saved.
%
% Outputs:
%   Writes a PNG file to disk and closes the figure.
    global params
    if ~exist(exportPath,'dir'), mkdir(exportPath); end
    fullP = fullfile(exportPath, [filenameStem '.png']);
    exportgraphics(figHandle, fullP, 'Resolution', params.exportRes);
    close(figHandle);
    fprintf('✔ Saved figure: %s\n', fullP);
end

function cmap = plasma(n)
% PLASMA Returns an interpolated plasma-style colormap.
%
% Inputs:
%   n - Number of colours requested.
%
% Outputs:
%   cmap - n-by-3 RGB colormap array.
    base = [ ...
        0.0504 0.0298 0.5280
        0.2546 0.0139 0.6154
        0.4176 0.0006 0.6584
        0.5627 0.0515 0.6415
        0.6928 0.1651 0.5645
        0.7982 0.2802 0.4695
        0.8814 0.3925 0.3832
        0.9440 0.5532 0.2871
        0.9796 0.7048 0.2129
        0.9892 0.8463 0.1400];
    x = linspace(0,1,size(base,1));
    xi = linspace(0,1,n);
    cmap = interp1(x, base, xi, 'linear');
end

function cmap = magma(n)
% MAGMA Returns an interpolated magma-style colormap.
%
% Inputs:
%   n - Number of colours requested. Defaults to 256 if omitted.
%
% Outputs:
%   cmap - n-by-3 RGB colormap array.
    if nargin < 1, n = 256; end
    base = [ ...
        0.0015 0.0005 0.0139
        0.0910 0.0277 0.2184
        0.2297 0.0599 0.4370
        0.3904 0.1004 0.5011
        0.5503 0.1612 0.5057
        0.7164 0.2140 0.4753
        0.8688 0.2877 0.4093
        0.9697 0.4761 0.2650
        0.9924 0.6800 0.2020
        0.9871 0.9914 0.7495];
    x = linspace(0,1,size(base,1));
    xi = linspace(0,1,n);
    cmap = interp1(x, base, xi, 'linear');
end

function pool = startParallelPool()
% STARTPARALLELPOOL Starts a clean local MATLAB parallel pool.
%
% Outputs:
%   pool - Active parallel pool object.
%
% Notes:
%   Existing pools are reused. Stale local cluster jobs and temporary job
%   storage are removed before creating a new pool.
    pool = gcp('nocreate');
    if ~isempty(pool)
        return;
    end
    myCluster = parcluster('local');
    jobStorage = fullfile(tempdir, 'MATLAB_Parallel_Reset');
    if exist(jobStorage, 'dir')
        try
            rmdir(jobStorage, 's');
        catch
        end
    end
    mkdir(jobStorage);
    myCluster.JobStorageLocation = jobStorage;
    jobs = myCluster.Jobs;
    if ~isempty(jobs)
        delete(jobs);
    end
    pool = parpool(myCluster);
end

function cleanupParallelPool()
% CLEANUPPARALLELPOOL Deletes active pools and clears stale parallel jobs.
%
% Outputs:
%   None.
%
% Notes:
%   This function is used before startup, after completion, and before fatal
%   errors to reduce problems caused by stale MATLAB parallel job storage.
    pool = gcp('nocreate');
    if ~isempty(pool)
        delete(pool);
    end
    myCluster = parcluster('local');
    jobs = myCluster.Jobs;
    if ~isempty(jobs)
        delete(jobs);
    end
    jobStorage = fullfile(tempdir, 'MATLAB_Parallel_Reset');
    if exist(jobStorage, 'dir')
        try
            rmdir(jobStorage, 's');
        catch
        end
    end
end
