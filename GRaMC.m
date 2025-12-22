%% GRaMC – Grain Reconstruction and Multi-stage Cleaning
% =========================================================================
% Description:
%   GRaMC is the companion script to MAPClean. It takes EBSD data (usually 
%   cleaned) and performs advanced grain reconstruction. It features:
%     1. Initial reconstruction based on misorientation angles.
%     2. Twin Boundary Merging (specifically for Anorthite laws).
%     3. Fake Inclusion Removal (merging small inclusions into hosts).
%     4. Background/Matrix Removal (isolating specific grain populations).
%
% Dependencies:
%   - MTEX Toolbox (v5.7.x or newer)
%
% Author: Rahul Subbaraman, University of Manchester
% Date: December 2025
% =========================================================================

clc; clear; close all;
import mtex.*;
setMTEXpref('generatingHelpMode','silent');
warning('off','all');

%% ================= USER CONFIGURATION =================

%% --- 1. Directories ---
dataDir       = fullfile(pwd, 'DataFiles');
exportDir     = fullfile(pwd, 'exports/GrainClean');
checkpointDir = fullfile(pwd, 'checkpoints');

% Create directories if they do not exist
if ~exist(exportDir,'dir'), mkdir(exportDir); end
if ~exist(checkpointDir,'dir'), mkdir(checkpointDir); end

%% --- 2. Stage Control Flags ---
% Set to 'false' to skip processing and load the relevant checkpoint instead.
run_initial     = false;   % Initial geometric reconstruction
run_twinMerge   = false;   % Merge twin boundaries (Anorthite)
run_fakeInc     = false;   % Merge small inclusions into host grains
run_bgRemoval   = false;   % Isolate large grains (Background selection)
run_Final       = true;    % Save final output

%% --- 3. Global Parameters ---
global params
params.exportRes        = 300;          % PNG export resolution (DPI)
params.grainThreshold   = 10*degree;    % Angle threshold for grain reconstruction
params.minPixelsGrains  = 5;            % Minimum pixel count to constitute a grain
params.bgSizeThres      = 300;          % Threshold to define "Background" grains (µm)
params.fakeIncSizeThres = 100;          % Threshold for "Fake Inclusion" removal (µm)

%% ================= MAIN PIPELINE =================

% File Selection (Targeting Anorthite files)
fileList = dir(fullfile(dataDir, '*_Anorthite.ctf'));
if isempty(fileList)
    error('GRaMC Error: No files matching "*_Anorthite.ctf" found in %s', dataDir);
end

for fi = 1:numel(fileList)
    % --- File Setup ---
    filePath = fullfile(dataDir, fileList(fi).name);
    [~, sampleName, ~] = fileparts(fileList(fi).name);
    sampleName = erase(sampleName, '_Anorthite'); % Clean naming
    
    % Setup export paths
    phaseExportPath = fullfile(exportDir, sampleName);
    if ~exist(phaseExportPath,'dir'), mkdir(phaseExportPath); end
    
    % Start Logging
    diaryFile = fullfile(phaseExportPath, [sampleName '_logfile.txt']);
    if exist(diaryFile, 'file'), delete(diaryFile); end
    diary(diaryFile); diary on;
    
    fprintf('\n===== Processing Sample: %s =====\n', sampleName);
    
    % --- Checkpoint Paths ---
    initialFile = fullfile(checkpointDir, sprintf('%s_initialGrains.mat', sampleName));
    twinFile    = fullfile(checkpointDir, sprintf('%s_TwinMergedGrains.mat', sampleName));
    incFile     = fullfile(checkpointDir, sprintf('%s_fakeIncRemGrains.mat', sampleName));
    bgFile      = fullfile(checkpointDir, sprintf('%s_bgRemGrains.mat', sampleName));
    FinalFile   = fullfile(checkpointDir, sprintf('%s_finalGrains.mat', sampleName));

    %% --- Load EBSD Data ---
    fprintf('Loading EBSD data...\n');
    ebsd = EBSD.load(filePath, 'convertSpatial2EulerReferenceFrame').gridify;
    fprintf('✔ EBSD loaded: %d points\n', numel(ebsd));
    ebsd_clean = ebsd;

    %% --- STEP 1: Initial Grain Reconstruction ---
    if run_initial
        fprintf('\n[STEP 1] Computing initial grains (Threshold: %.1f°)...\n', params.grainThreshold/degree);
        [ogGrains, ebsd_clean.grainId] = calcGrains(ebsd_clean, 'angle', params.grainThreshold);
        
        % Filter small noise grains
        smallMask = ogGrains.grainSize < params.minPixelsGrains;
        if any(smallMask)
            fprintf('  -> Removing %d small noise grains (<%d pixels).\n', sum(smallMask), params.minPixelsGrains);
            ogGrains(smallMask) = []; 
            % Update EBSD grid to reflect removal
            ebsd_clean(ismember(ebsd_clean.grainId, find(smallMask))) = []; 
            ebsd_clean = ebsd_clean.gridify;
        end
        
        save(initialFile, 'ogGrains', 'filePath');
        fprintf('✔ Checkpoint saved: %s\n', initialFile);
        plotGrainMaps(ogGrains, phaseExportPath, sampleName, '01_PreMerge');
    else
        checkExists(initialFile, 'STEP 1');
        S = load(initialFile, 'ogGrains', 'filePath');
        ogGrains = S.ogGrains;
        fprintf('✔ [STEP 1] Loaded checkpoint: %d grains.\n', ogGrains.length);
    end

    %% --- STEP 2: Twin Boundary Merging ---
    if run_twinMerge
        fprintf('\n[STEP 2] Computing Twin Merges (Anorthite Laws)...\n');
        gB = ogGrains.boundary;
        
        % Define Twin Laws
        twinLaws = {
            {'Albite',    orientation.byAxisAngle(vector3d(0,1,0), 180*degree, ebsd_clean.CS), 5*degree},
            {'Pericline', orientation.byAxisAngle(vector3d(1,0,0), 180*degree, ebsd_clean.CS), 5*degree},
            {'Carlsbad',  orientation.byAxisAngle(vector3d(0,0,1), 180*degree, ebsd_clean.CS), 5*degree},
            {'Manebach',  orientation(reflection(Miller(0,0,1,ebsd_clean.CS))), 5*degree},
            {'Baveno',    orientation(reflection(Miller(0,2,1,ebsd_clean.CS))), 5*degree} 
        };
        
        twinId = zeros(size(gB));
        for i = 1:length(twinLaws)
            law = twinLaws{i};
            isTwin = angle(gB.misorientation, law{2}) < law{3};
            newTwins = isTwin & (twinId == 0);
            twinId(newTwins) = i;
            fprintf('  -> Detected %-10s boundaries: %d\n', law{1}, sum(newTwins));
        end
        
        % Merge Grains across Twin Boundaries
        twinBoundaries = gB(twinId > 0);
        gidpair = unique(sort(twinBoundaries.grainId, 2), 'rows');
        [TwinMergedGrains, parentId] = mergeGrainsAndComputeGOS(ogGrains, gidpair);
        
        save(twinFile, 'TwinMergedGrains', 'parentId', 'filePath');
        fprintf('✔ Checkpoint saved: %s\n', twinFile);
        plotGrainMaps(TwinMergedGrains, phaseExportPath, sampleName, '02_TwinMerge');
    else
        checkExists(twinFile, 'STEP 2');
        S = load(twinFile, 'TwinMergedGrains', 'parentId', 'filePath');
        TwinMergedGrains = S.TwinMergedGrains;
        parentId = S.parentId;
        fprintf('✔ [STEP 2] Loaded checkpoint: %d grains.\n', TwinMergedGrains.length);
    end

    %% --- STEP 3: Fake Inclusion Removal ---
    if run_fakeInc 
        fprintf('\n[STEP 3] Removing Fake Inclusions (Threshold: %d µm)...\n', params.fakeIncSizeThres);
        stepSize = (max(ebsd_clean.x(:)) - min(ebsd_clean.x(:))) / (ebsd.size(2) - 1);
        
        [isIncl, hostId] = TwinMergedGrains.isInclusion;
        gosMask  = TwinMergedGrains.GOS == 0; 
        % Convert micron threshold to pixels
        sizeMask = TwinMergedGrains.grainSize <= ceil(params.fakeIncSizeThres / stepSize);
        
        % Identify fake inclusions (Small/Zero-GOS grains trapped inside others)
        fakeIncMask = (gosMask | sizeMask) & isIncl;
        
        if any(fakeIncMask)
            fprintf('  -> Merging %d fake inclusions into host grains.\n', sum(fakeIncMask));
            gidpair = [find(fakeIncMask) hostId(fakeIncMask)];
            [fakeIncRemGrains, parentMap] = mergeGrainsAndComputeGOS(TwinMergedGrains, gidpair);
        else
            fprintf('  -> No fake inclusions detected.\n');
            fakeIncRemGrains = TwinMergedGrains;
            parentMap = [];
        end
        
        save(incFile, 'fakeIncRemGrains', 'parentMap', 'filePath');
        fprintf('✔ Checkpoint saved: %s\n', incFile);
        plotGrainMaps(fakeIncRemGrains, phaseExportPath, sampleName, '03_FakeIncFix');
    else
        checkExists(incFile, 'STEP 3');
        S = load(incFile, 'fakeIncRemGrains', 'parentMap', 'filePath');
        fakeIncRemGrains = S.fakeIncRemGrains;
        parentMap = S.parentMap;
        fprintf('✔ [STEP 3] Loaded checkpoint: %d grains.\n', fakeIncRemGrains.length);
    end

    %% --- STEP 4: Background Selection (Matrix Removal) ---
    if run_bgRemoval
        fprintf('\n[STEP 4] Isolating Large Grains (Background Selection > %d µm)...\n', params.bgSizeThres);
        [isIncl, ~] = fakeIncRemGrains.isInclusion;
        stepSize = (max(ebsd_clean.x(:)) - min(ebsd_clean.x(:))) / (ebsd.size(2) - 1);
        
        % Select Large Grains (The "Background" or Main Phase)
        sizeMask = fakeIncRemGrains.grainSize > ceil(params.bgSizeThres / stepSize);
        bgMask = sizeMask & ~isIncl; % Must be large and not an inclusion
        
        % NOTE: 'bgRemGrains' here refers to the retained large grains
        bgRemGrains = fakeIncRemGrains(bgMask);
        fprintf('  -> Retained %d large grains.\n', bgRemGrains.length);
        
        save(bgFile, 'bgRemGrains', 'filePath');
        fprintf('✔ Checkpoint saved: %s\n', bgFile);
        plotGrainMaps(bgRemGrains, phaseExportPath, sampleName, '04_FinalSelection');
    else
        checkExists(bgFile, 'STEP 4');
        S = load(bgFile, 'bgRemGrains', 'filePath');
        bgRemGrains = S.bgRemGrains;
        fprintf('✔ [STEP 4] Loaded checkpoint: %d grains.\n', bgRemGrains.length);
    end

    %% --- STEP 5: Final Output ---
    if run_Final 
        finalGrains = bgRemGrains;
        save(FinalFile, 'finalGrains', 'filePath');
        fprintf('\n[FINAL] Grains saved to: %s\n', FinalFile);
    end
    
    diary off;
end

%% ================= HELPER FUNCTIONS =================

function checkExists(fileP, stepName)
    if ~exist(fileP, 'file')
        error('%s disabled but checkpoint file missing: %s', stepName, fileP);
    end
end

function [mergedGrains, parentId] = mergeGrainsAndComputeGOS(grains, gidpair)
% MERGEGRAINSANDCOMPUTEGOS Merges grains and recalculates properties.
    % 1. Merge grains
    [mergedGrains, parentId] = merge(grains, gidpair);
    
    % 2. Parent Orientation: Assign orientation of largest child grain
    [~, sortIdx] = sort(grains.area, 'descend');
    sortedParentId = parentId(sortIdx);
    [uParents, uIdx] = unique(sortedParentId, 'stable');
    bestChildIdx = sortIdx(uIdx);
    
    newOri = mergedGrains.meanOrientation;  
    newOri(uParents) = grains.meanOrientation(bestChildIdx);
    mergedGrains.meanOrientation = newOri;
    
    % 3. Parent GOS: Area-weighted average
    weightedSum = accumarray(parentId, grains.GOS .* grains.area);
    totalArea   = accumarray(parentId, grains.area);
    mergedGrains.prop.GOS = weightedSum ./ totalArea;
end

function plotGrainMaps(grains, exportPath, sampleName, keyWord)    
    global params
    % 1. Boundaries
    fBound = figure('Visible','off');
    plot(grains.boundary, 'lineColor','k'); axis equal tight;
    savePNG(fBound, sprintf('%s_%s_Boundaries', sampleName, keyWord), exportPath);
    
    % 2. GOS
    fGOS = figure('Visible','off');
    plot(grains, grains.prop.GOS./degree, 'backgroundColor','k');
    mtexColorbar; colormap(gca, flipud(magma));
    setColorRange([0 max(grains.prop.GOS./degree)]); axis equal tight;
    savePNG(fGOS, sprintf('%s_%s_GOS', sampleName, keyWord), exportPath);
    
    % 3. IPF (Z-direction)
    ipfKey = ipfColorKey(grains);
    ipfKey.inversePoleFigureDirection = vector3d.Z;
    colors = ipfKey.orientation2color(grains.meanOrientation);
    fIPF = figure('Visible','off');
    plot(grains, colors);
    hold on; plot(grains.boundary, 'lineColor','k', 'lineWidth',1); hold off;
    axis equal tight;
    savePNG(fIPF, sprintf('%s_%s_IPF', sampleName, keyWord), exportPath);
end

function savePNG(figHandle, filenameStem, exportPath)
    global params
    if ~exist(exportPath,'dir'), mkdir(exportPath); end
    fullP = fullfile(exportPath, [filenameStem '.png']);
    exportgraphics(figHandle, fullP, 'Resolution', params.exportRes);
    close(figHandle);
end
