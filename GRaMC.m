%% GRaMC – Grain Reconstruction and Multi-stage Cleaning
% =========================================================================
% Description:
%   GRaMC is a grain-scale downstream pipeline for MAPClean outputs.
%   It reconstructs grains from phase-specific EBSD data, removes
%   topologically unreliable GOS = 0 grains, merges crystallographically
%   compatible twin-related daughter grains, and optionally applies a final
%   area-based population filter before grain statistics.
%
%   This version performs:
%     1. Phase splitting from MAPClean outputs
%     2. Initial grain reconstruction on phase-specific Anorthite data
%     3. Removal of GOS = 0 grains by resetting their pixels to notIndexed
%     4. Twin boundary detection and twin merging
%     5. Pixel-level EBSD reorientation into the merged parent frame
%     6. Grain-size histogram visualisation in log10(area) space
%     7. Area-based filtering of unwanted small-grain populations
%
%   Outputs:
%     - Boundary, GOS, and IPF maps at each major stage
%     - Extra twin-boundary classification map at twin stage
%     - Extra area-distribution histogram and colour-linked grain map
%     - Stage-wise .mat checkpoints containing both EBSD and grains
%
% Dependencies:
%   - MTEX Toolbox (tested on v 6.0.0)
%
% Author: Rahul Subbaraman
% Version 1.0: December 2025, 2.0: April 2026
% =========================================================================

clc; clear; close all;
import mtex.*;
setMTEXpref('generatingHelpMode','silent');
warning('off','all');

%% ================= USER CONFIGURATION =================

%% --- 1. Directories ---
dataDir       = fullfile(pwd, 'DataFiles');
exportDir     = fullfile(pwd, 'GRaMC');
checkpointDir = fullfile(pwd, 'checkpoints');

if ~exist(exportDir,'dir'), mkdir(exportDir); end
if ~exist(checkpointDir,'dir'), mkdir(checkpointDir); end

%% --- 2. Stage Control Flags ---
runPhaseSplit  = true;   % Split MAPClean output into phase-specific .ctf files
run_initial    = true;   % Initial geometric reconstruction
run_gosClean   = true;   % Remove GOS = 0 grains
run_twinMerge  = true;   % Merge twin boundaries and update EBSD pixels
run_sizeThresh = true;   % Create colour-coded grain size distribution
run_areaFilter = true;    % Area filter 
run_Final      = true;   % Save final output

%% --- 3. Global Parameters ---
global params
params.exportRes       = 300;          % PNG export resolution (DPI)
params.grainThreshold  = 10*degree;    % Angle threshold for grain reconstruction
params.minPixelsGrains = 5;            % Minimum pixel count to constitute a grain
params.binsPerDecade   = 10;            % Hitogram bins per 10^x to 10^(x+1)
params.logAreaThresh = containers.Map('KeyType','char','ValueType','double'); % Manual log10(area) threshold per sample
params.logAreaThresh('01a6') = 4.4;   % set manually from histogram
params.logAreaThresh('4N3C') = 4.4;   % set manually from histogram
params.logAreaThresh('024')  = 4.4;   % set manually from histogram
params.logAreaThresh('01a2') = 4.6;   % set manually from histogram

%% ================= MAIN PIPELINE =================
sampleNames = {'01a6','4N3C','024','01a2'}; % samplename must match "filename_clean.ctf" without "_clean.ctf"
fileList = cellfun(@(s) dir(fullfile(dataDir,[s '_clean.ctf'])), sampleNames, 'UniformOutput', false);
fileList = vertcat(fileList{:});

if isempty(fileList)
    error('GRaMC Error: No *_clean.ctf files found in %s. Run MAPClean first.', dataDir);
end

for fi = 1:numel(fileList)
    tSample = tic;

    %% --- File Setup ---
    filePath = fullfile(dataDir, fileList(fi).name);
    [~, sampleName, ~] = fileparts(fileList(fi).name);
    sampleName = erase(sampleName, '_clean');
    sampleName = erase(sampleName, '_Anorthite');

    phaseExportPath = fullfile(exportDir, sampleName);
    if ~exist(phaseExportPath,'dir'), mkdir(phaseExportPath); end

    %% --- Logging ---
    diaryFile = fullfile(phaseExportPath, [sampleName '_logfile.txt']);
    diary off;
    if exist(diaryFile, 'file'), delete(diaryFile); end
    diary(diaryFile); diary on;

    fprintf('\n===== Processing Sample: %s =====\n', sampleName);

    %% --- Checkpoint Paths ---
    splitFile   = fullfile(checkpointDir, sprintf('%s_phaseSplit.mat', sampleName));
    initialFile = fullfile(checkpointDir, sprintf('%s_initialGrains.mat', sampleName));
    gosFile     = fullfile(checkpointDir, sprintf('%s_gosCleanGrains.mat', sampleName));
    twinFile    = fullfile(checkpointDir, sprintf('%s_TwinMergedGrains.mat', sampleName));
    areaFile = fullfile(checkpointDir, sprintf('%s_areaFilteredGrains.mat', sampleName));
    FinalFile   = fullfile(checkpointDir, sprintf('%s_finalGrains.mat', sampleName));

    %% --- Load Full MAPClean Output ---
    tLoad = tic;
    fprintf('Loading full cleaned EBSD data...\n');
    ebsd_full = EBSD.load(filePath, 'convertSpatial2EulerReferenceFrame').gridify;
    fprintf('✔ Full EBSD loaded: %d points (%.2f s)\n', numel(ebsd_full), toc(tLoad));

    %% --- STEP 0: Phase Split Export / Check ---
    % Splits the full MAPClean output into phase-specific .ctf files so
    % downstream reconstruction can be performed on one mineral population
    % at a time. If phase splitting is disabled, the code instead checks
    % that the required phase-specific files already exist before continuing.
    tStep0 = tic;
    fprintf('\n[STEP 0] Phase Split Export / Check...\n');

    phases = ebsd_full.mineralList;
    notIndexedId = find(strcmpi(phases,'notIndexed'));
    MinPhaseIds = setdiff(1:numel(phases), notIndexedId);
    MinPhaseNames = phases(MinPhaseIds);
    [Nrow, Ncol] = ebsd_full.size;

    if runPhaseSplit
        for p = 1:numel(MinPhaseNames)
            phaseName = MinPhaseNames{p};
            phaseIdx = find(strcmp(phases, phaseName));

            ebsd_phase = ebsd_full;
            maskLinear = ebsd_phase.phaseId ~= phaseIdx;
            maskGrid = reshape(maskLinear, Nrow, Ncol);
            ebsd_phase(maskGrid).phaseId = notIndexedId;

            validCount = sum(ebsd_phase.phaseId == phaseIdx);
            fprintf('  -> %s has %d indexed points.\n', phaseName, validCount);

            outFile = fullfile(dataDir, sprintf('%s_%s.ctf', sampleName, phaseName));
            ebsd_phase.export(outFile);
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
            error('Phase files do not exist. Set runPhaseSplit = true and rerun.');
        end

        fprintf('✔ [STEP 0] All phase files exist. Continuing (%.2f s)\n', toc(tStep0));
    end

    %% --- Load Phase-Specific EBSD for GRaMC (Anorthite)---
    tPhaseLoad = tic;
    phaseFile = fullfile(dataDir, sprintf('%s_Anorthite.ctf', sampleName));
    fprintf('Loading phase-specific EBSD: %s\n', phaseFile);
    ebsd_phase = EBSD.load(phaseFile, 'convertSpatial2EulerReferenceFrame').gridify;
    fprintf('✔ Anorthite EBSD loaded: %d points (%.2f s)\n', numel(ebsd_phase), toc(tPhaseLoad));

    phases_phase = ebsd_phase.mineralList;
    notIndexedId_phase = find(strcmpi(phases_phase,'notIndexed'));
    anorthiteId = find(strcmpi(phases_phase,'Anorthite'));
    if isempty(anorthiteId)
        error('Anorthite phase not found in %s', phaseFile);
    end

    %% --- STEP 1: Initial Grain Reconstruction ---
    % Reconstructs grains from the phase-specific EBSD map using the global
    % grain-threshold angle. A minimum pixel-count filter is then applied
    % to remove sub-resolution noise grains. Pixels belonging to removed
    % grains are reset to notIndexed, after which grains are recalculated
    % so that the EBSD object and grain object remain fully consistent.
    tStep1 = tic;
    if run_initial
        fprintf('\n[STEP 1] Computing initial grains (Threshold: %.1f°)...\n', params.grainThreshold/degree);

        ebsd_initial = ebsd_phase;
        [grains_initial, ebsd_initial.grainId] = calcGrains(ebsd_initial, 'angle', params.grainThreshold);
        fprintf('  -> Initial grains found: %d\n', grains_initial.length);

        smallMask = grains_initial.grainSize < params.minPixelsGrains;
        if any(smallMask)
            fprintf('  -> Removing %d small noise grains (<%d pixels).\n', sum(smallMask), params.minPixelsGrains);

            grainIdsToRemove = grains_initial.id(smallMask);
            pixelMaskRemove = ismember(ebsd_initial.grainId, grainIdsToRemove);
            ebsd_initial(pixelMaskRemove).phaseId = notIndexedId_phase;

            [grains_initial, ebsd_initial.grainId] = calcGrains(ebsd_initial, 'angle', params.grainThreshold);
            fprintf('  -> Grain count after small-grain removal: %d\n', grains_initial.length);
        end

        save(initialFile, 'ebsd_initial', 'grains_initial', 'filePath');
        fprintf('✔ Checkpoint saved: %s\n', initialFile);

        plotGrainMaps(grains_initial, phaseExportPath, sampleName, '01_Initial');
        fprintf('✔ [STEP 1] Complete (%.2f s)\n', toc(tStep1));
    else
        checkExists(initialFile, 'STEP 1');
        S = load(initialFile, 'ebsd_initial', 'grains_initial', 'filePath');
        ebsd_initial = S.ebsd_initial;
        grains_initial = S.grains_initial;
        fprintf('✔ [STEP 1] Initial grain reconstruction checkpoint loaded: %d grains, in %.2f s\n', grains_initial.length, toc(tStep1));
    end

    %% --- STEP 2: GOS = 0 Grain Removal ---
    % Identifies grains with zero Grain Orientation Spread (GOS), which in
    % this workflow are treated as non-physical or topologically unreliable
    % grain objects. Pixels belonging to these grains are reset to
    % notIndexed and the grain structure is recalculated. A diagnostic map
    % is also exported to show the spatial distribution of GOS = 0 grains.
    tStep2 = tic;
    if run_gosClean
        fprintf('\n[STEP 2] Removing GOS = 0 grains...\n');

        ebsd_gos = ebsd_initial;
        grains_gos_pre = grains_initial;

        gosMask = grains_gos_pre.GOS == 0;
        nGOS0 = sum(gosMask);
        fprintf('  -> GOS = 0 grains found: %d\n', nGOS0);

        plotGOSZeroGrains(grains_gos_pre, gosMask, phaseExportPath, sampleName, '02a_GOS0_Grains');

        if nGOS0 > 0
            grainIdsToRemove = grains_gos_pre.id(gosMask);
            pixelMaskRemove = ismember(ebsd_gos.grainId, grainIdsToRemove);

            ebsd_gos(pixelMaskRemove).phaseId = notIndexedId_phase;
            [grains_gosclean, ebsd_gos.grainId] = calcGrains(ebsd_gos, 'angle', params.grainThreshold);
        else
            grains_gosclean = grains_gos_pre;
        end

        fprintf('  -> Grain count after GOS = 0 removal: %d\n', grains_gosclean.length);

        save(gosFile, 'ebsd_gos', 'grains_gosclean', 'gosMask', 'filePath');
        fprintf('✔ Checkpoint saved: %s\n', gosFile);

        plotGrainMaps(grains_gosclean, phaseExportPath, sampleName, '02_GOSClean');
        fprintf('✔ [STEP 2] Complete (%.2f s)\n', toc(tStep2));
    else
        checkExists(gosFile, 'STEP 2');
        S = load(gosFile, 'ebsd_gos', 'grains_gosclean', 'gosMask', 'filePath');
        ebsd_gos = S.ebsd_gos;
        grains_gosclean = S.grains_gosclean;
        gosMask = S.gosMask;
        fprintf('✔ [STEP 2] Grains with GOS = 0 removal checkpoint loaded: %d grains, in %.2f s\n', grains_gosclean.length, toc(tStep2));
    end

    %% --- STEP 3: Twin Boundary Merging + EBSD Update ---
    % Detects crystallographically valid Anorthite twin boundaries,
    % classifies them by twin law, and merges daughter grains linked by
    % those boundaries into parent grains. The largest-area daughter is
    % retained as the reference orientation for each merged parent.
    %
    % Pixel-level EBSD orientations are then updated so that all rotating
    % daughter grains are brought into the parent reference frame. This is
    % done parent-by-parent, with only non-identity daughter grains being
    % rotated. All daughter pixels are reassigned to the merged parent
    % grain ID, after which grains are recalculated from the updated EBSD.
    tStep3 = tic;
    if run_twinMerge
        fprintf('\n[STEP 3] Computing Twin Merges (Anorthite Laws)...\n');

        ebsd_twin = ebsd_gos;
        grains_preTwin = grains_gosclean;
        nGrainsBefore = grains_preTwin.length;
        fprintf('  -> Grain count before twin merge: %d\n', nGrainsBefore);

        gB = grains_preTwin.boundary;
        cs = ebsd_twin('Anorthite').CS;

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
        plotTwinBoundaryMap(grains_preTwin, gB, twinId, twinLaws, phaseExportPath, sampleName, '03a_TwinClassification');
        twinBoundaries = gB(twinId > 0);
        gidpair = unique(sort(twinBoundaries.grainId, 2), 'rows');

        if isempty(gidpair)
            fprintf('  -> No twin-linked grain pairs found. Skipping merge.\n');
            grains_twin = grains_preTwin;
            fprintf('  -> Grain count after twin merge: %d\n', grains_twin.length);
            save(twinFile, 'ebsd_twin', 'grains_twin', 'twinId', 'filePath');
            fprintf('✔ Checkpoint saved: %s\n', twinFile);
            plotGrainMaps(grains_twin, phaseExportPath, sampleName, '03_TwinMerge');
            fprintf('✔ [STEP 3] Complete (%.2f s)\n', toc(tStep3));
        else
            tMerge = tic;
            [grains_twin_raw, parentId] = mergeGrainsAndComputeGOS(grains_preTwin, gidpair);
            fprintf('  -> Grain-object merge complete (%.2f s)\n', toc(tMerge));
            % Precompute pixel bookkeeping
            oldGrainId = ebsd_twin.grainId;
            validPix = oldGrainId > 0 & oldGrainId <= numel(parentId);
            validPixIdx = find(validPix);
            pixByPreId = accumarray(oldGrainId(validPix), validPixIdx, [numel(parentId), 1], @(x){x});
            % Unique post IDs and reference daughters 
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
            % Parallel pool
            pool = gcp('nocreate');
            if isempty(pool)
                pool = parpool('local');
            end
            fprintf('  -> Parallel pool active with %d workers\n', pool.NumWorkers);
            % Process parent grains in chunks
            rotatedPixelTotal = 0;
            preMeanOri = grains_preTwin.meanOrientation;
            postMeanOri = grains_twin_raw.meanOrientation;
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
                    parfor kk = 1:nThis
                        thisPostId = thisPostChunk(kk);
                        daughters = find(postIdPerPreId == thisPostId);
                        refDaughter = refDaughterByPost(thisPostId);
                        rotDaughters = daughters(daughters ~= refDaughter);
                        pixIdxAll = [];
                        quatAll = zeros(0,4);
                        postIdAll = [];
                        qRef = quaternion(preMeanOri(refDaughter));
                        for dd = 1:numel(rotDaughters)
                            gid = rotDaughters(dd);
                            pixIdx = pixByPreId{gid};
                            if isempty(pixIdx)
                                continue;
                            end
                            qD = quaternion(preMeanOri(gid));
                            rotOp = qRef * inv(qD);
                            qPix = quaternion(ebsd_twin(pixIdx).orientations);
                            qRot = rotOp * qPix;
                            pixIdxAll = [pixIdxAll; pixIdx(:)]; %#ok<AGROW>
                            quatAll = [quatAll; [qRot.a(:), qRot.b(:), qRot.c(:), qRot.d(:)]]; %#ok<AGROW>
                            postIdAll = [postIdAll; repmat(thisPostId, numel(pixIdx), 1)]; %#ok<AGROW>
                        end
                        % also remap reference daughter pixels to post ID, no rotation
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
                    % Serial write-back
                    for kk = 1:nThis
                        pixIdxAll = chunkPixIdx{kk};
                        quatAll   = chunkQuat{kk};
                        postIdAll = chunkPostId{kk};
                        if isempty(pixIdxAll)
                            continue;
                        end
                        % rotate only rows with nonzero quaternion update
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
            % Recalculate grains from updated EBSD
            tRecalc = tic;
            [grains_twin, ebsd_twin.grainId] = calcGrains(ebsd_twin, 'angle', params.grainThreshold);
            fprintf('  -> Recalculated grains from updated EBSD (%.2f s)\n', toc(tRecalc));
            fprintf('  -> Grain count after twin merge: %d\n', grains_twin.length);
            save(twinFile, 'ebsd_twin', 'grains_twin', 'twinId', 'filePath');
            fprintf('✔ Checkpoint saved: %s\n', twinFile);
            plotGrainMaps(grains_twin, phaseExportPath, sampleName, '03_TwinMerge');
            fprintf('✔ [STEP 3] Complete (%.2f s)\n', toc(tStep3));
        end
    else
        checkExists(twinFile, 'STEP 3');
        S = load(twinFile, 'ebsd_twin', 'grains_twin', 'twinId', 'filePath');
        ebsd_twin = S.ebsd_twin;
        grains_twin = S.grains_twin;
        twinId = S.twinId;
        fprintf('✔ [STEP 3] Twin merge checkpoint loaded: %d grains, in %.2f s\n', grains_twin.length, toc(tStep3));
    end

    %% --- STEP 4: log10(Area) Histogram + Matching Grain-Bin Map ---
    % Builds a diagnostic grain-size distribution in log10(area) space and
    % exports a matching spatial grain map using the same colour bins. This
    % step is intended purely for visual threshold selection: it allows the
    % user to correlate histogram bins with actual grain populations in the
    % map before applying any area-based filtering.
    tStep4 = tic;
    if run_sizeThresh
        fprintf('\n[STEP 4] log10(area) histogram + matching grain-bin map...\n');
        areas = grains_twin.area;
        logAreas = log10(areas);
        % --- User control: bins per decade ---
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
        % Histogram
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
        title(sprintf('%s: log10(Grain Area) Distribution', sampleName), 'Interpreter', 'none');
        grid on;
        hold off;
        savePNG(fHist, sprintf('%s_04_LogAreaHistogram', sampleName), phaseExportPath);

        % Grain map with same bin colours
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
        title(sprintf('%s: Grain Area Bin Map', sampleName), 'Interpreter', 'none');
        hold off;

        colormap(cmap);
        cb = colorbar;
        cb.Ticks = linspace(1/(2*nbins), 1 - 1/(2*nbins), nbins);
        cb.TickLabels = arrayfun(@(i) sprintf('%.2f–%.2f', edges(i), edges(i+1)), 1:nbins, 'UniformOutput', false);
        cb.Label.String = 'log10(grain area) bins';
        cb.FontSize = 6;
        cb.Label.FontSize = 6;

        savePNG(fMap, sprintf('%s_04_LogAreaBinMap', sampleName), phaseExportPath);

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
    % Applies a user-defined log10(area) threshold to remove unwanted
    % small-grain populations, such as microlites or other grains outside
    % the target size range for final grain statistics. Grains below the
    % threshold are removed by resetting their pixels to notIndexed, and
    % grains are then recalculated so the EBSD and grain objects remain
    % consistent. If the threshold is set to zero for a sample, this step
    % is skipped and the pre-filtered reconstruction is retained.
    tStep5 = tic;
    if run_areaFilter
        if ~isKey(params.logAreaThresh, sampleName)
            error('No log10(area) threshold defined for sample %s', sampleName);
        end

        logAreaThresh = params.logAreaThresh(sampleName);

        if logAreaThresh == 0
            fprintf('\n[STEP 5] Area filter skipped because log10(area) threshold = 0 for %s\n', sampleName);
            ebsd_area = ebsd_twin;
            grains_area = grains_twin;
            fprintf('✔ [STEP 5] Skipped (%.2f s)\n', toc(tStep5));
        else
            fprintf('\n[STEP 5] Area-based grain filtering...\n');

            areaThresh = 10^logAreaThresh;

            fprintf('  -> log10(area) threshold for %s = %.6f\n', sampleName, logAreaThresh);
            fprintf('  -> area threshold for %s = %.6f\n', sampleName, areaThresh);

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

            save(areaFile, 'ebsd_area', 'grains_area', 'logAreaThresh', 'areaThresh', 'filePath');
            fprintf('✔ Checkpoint saved: %s\n', areaFile);

            plotGrainMaps(grains_area, phaseExportPath, sampleName, '05_AreaFiltered');

            fprintf('✔ [STEP 5] Complete (%.2f s)\n', toc(tStep5));
        end
    else
        checkExists(areaFile, 'STEP 5');
        S = load(areaFile, 'ebsd_area', 'grains_area', 'logAreaThresh', 'areaThresh', 'filePath');
        ebsd_area = S.ebsd_area;
        grains_area = S.grains_area;
        logAreaThresh = S.logAreaThresh;
        areaThresh = S.areaThresh;
        fprintf('✔ [STEP 5] Loaded checkpoint: %d grains.\n', grains_area.length);
        fprintf('✔ [STEP 5] Complete (%.2f s)\n', toc(tStep5));
    end

    %% --- STEP 6: Final Output ---
    % Selects the final EBSD/grain state based on the highest completed
    % processing stage and saves it as the final checkpoint for downstream
    % statistics, inspection, and export.
    tStep6 = tic;
    if run_Final
        if run_areaFilter && isKey(params.logAreaThresh, sampleName) && params.logAreaThresh(sampleName) ~= 0
            ebsd_final = ebsd_area;
            finalGrains = grains_area;
        elseif run_twinMerge || exist(twinFile, 'file')
            ebsd_final = ebsd_twin;
            finalGrains = grains_twin;
        elseif run_gosClean || exist(gosFile, 'file')
            ebsd_final = ebsd_gos;
            finalGrains = grains_gosclean;
        else
            ebsd_final = ebsd_initial;
            finalGrains = grains_initial;
        end

        save(FinalFile, 'ebsd_final', 'finalGrains', 'filePath');
        fprintf('\n[FINAL] Grains saved to: %s\n', FinalFile);
        fprintf('✔ [STEP 6] Complete (%.2f s)\n', toc(tStep6));
    else
        fprintf('✔ [STEP 6] Skipped (%.2f s)\n', toc(tStep6));
    end
end

%% ================= HELPER FUNCTIONS =================

function checkExists(fileP, stepName)
% CHECKEXISTS Throws an error if a required checkpoint file is missing.
    if ~exist(fileP, 'file')
        error('%s disabled but checkpoint file missing: %s', stepName, fileP);
    end
end

function [mergedGrains, parentId] = mergeGrainsAndComputeGOS(grains, gidpair)
% MERGEGRAINSANDCOMPUTEGOS Merges grains and recalculates parent properties.
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

function plotGrainMaps(grains, exportPath, sampleName, keyWord)
% PLOTGRAINMAPS Exports standard boundary, GOS, and IPF grain maps.
    global params
    % Grain Boundary
    fBound = figure('Visible','off');
    plot(grains.boundary, 'lineColor','k');
    axis equal tight;
    savePNG(fBound, sprintf('%s_%s_Boundaries', sampleName, keyWord), exportPath);
    % GOS
    fGOS = figure('Visible','off');
    plot(grains, grains.prop.GOS./degree, 'backgroundColor','k');
    mtexColorbar;
    colormap(gca, flipud(magma));
    setColorRange([0 max(grains.prop.GOS./degree)]);
    axis equal tight;
    savePNG(fGOS, sprintf('%s_%s_GOS', sampleName, keyWord), exportPath);
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
    savePNG(fIPF, sprintf('%s_%s_IPF', sampleName, keyWord), exportPath);
end

function plotGOSZeroGrains(grains, gosMask, exportPath, sampleName, keyWord)
% PLOTGOSZEROGRAINS Exports a diagnostic map of grains with GOS = 0.
    fGOS0 = figure('Visible','off');
    hold on;
    if any(gosMask)
        plot(grains(gosMask), 'FaceColor', [1 0 0], 'EdgeColor', 'none');
    end
    plot(grains.boundary, 'lineColor', 'k', 'lineWidth', 0.5);
    axis equal tight;
    hold off;
    savePNG(fGOS0, sprintf('%s_%s', sampleName, keyWord), exportPath);
end

function plotTwinBoundaryMap(~, gB, twinId, twinLaws, exportPath, sampleName, keyWord)
% PLOTTWINBOUNDARYMAP Exports a classified twin-boundary map by twin law.
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
    savePNG(fTwins, sprintf('%s_%s_TwinBoundaries', sampleName, keyWord), exportPath);
end

function savePNG(figHandle, filenameStem, exportPath)
% SAVEPNG Saves a figure to PNG at the global export resolution.
    global params
    if ~exist(exportPath,'dir'), mkdir(exportPath); end
    fullP = fullfile(exportPath, [filenameStem '.png']);
    exportgraphics(figHandle, fullP, 'Resolution', params.exportRes);
    close(figHandle);
    fprintf('✔ Saved figure: %s\n', fullP);
end

function cmap = plasma(n)
% PLASMA Returns an interpolated plasma-style colormap with n colours.
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

