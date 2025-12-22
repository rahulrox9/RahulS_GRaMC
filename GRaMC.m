%% GRaMC – Grain Reconstruction and Multi-stage Cleaning
% =========================================================================
% Description:
%   GRaMC is an advanced microstructural restoration engine designed to work downstream of MAPClean. While standard reconstruction simply groups pixels by misorientation, GRaMC applies crystallographic and topological heuristics to repair "over-segmented" datasets (e.g., Plagioclase).
%
%   Key Capabilities:
%     1. Crystallographic Repair: Automatically detects and heals twin 
%        boundaries (Albite, Pericline, etc.) to restore parent grains.
%     2. Topological Cleaning: Identifies "fake inclusions" using multi-
%        criteria logic (Size + GOS + Enclosure) and absorbs them into hosts.
%     3. Data Fidelity: Preserves statistical validity during merging by 
%        recalculating area-weighted GOS and dominant orientations.
%     4. Population Isolation: Filters groundmass to isolate phenocrysts.
%
% Dependencies:
%   - MTEX Toolbox (tested on v 6.0.0)
%
% Author: Rahul Subbaraman 
% Date: December 2025
% License: GPLv3
% =========================================================================

clc; clear; close all;
import mtex.*;
setMTEXpref('generatingHelpMode','silent');
warning('off','all');

%% ================= USER CONFIGURATION =================

%% --- 1. Directories ---
% Assumes "Unified Workspace" structure
dataDir       = fullfile(pwd, 'DataFiles');
exportDir     = fullfile(pwd, 'exports/GrainClean');
checkpointDir = fullfile(pwd, 'checkpoints');

% Create directories if they do not exist
if ~exist(exportDir,'dir'), mkdir(exportDir); end
if ~exist(checkpointDir,'dir'), mkdir(checkpointDir); end

%% --- 2. Stage Control Flags ---
% Set to 'false' to skip processing and load the relevant checkpoint instead.
run_initial     = false;   % Initial geometric reconstruction
run_twinMerge   = true;    % Merge twin boundaries (Anorthite)
run_fakeInc     = true;    % Merge small inclusions into host grains
run_bgRemoval   = true;    % Isolate large grains (Background selection)
run_Final       = true;    % Save final output

%% --- 3. Global Parameters ---
global params
params.exportRes        = 300;          % PNG export resolution (DPI)
params.grainThreshold   = 10*degree;    % Angle threshold for grain reconstruction
params.minPixelsGrains  = 5;            % Minimum pixel count to constitute a grain
params.bgSizeThres      = 300;          % Threshold to define "Background" grains (µm)
params.fakeIncSizeThres = 100;          % Threshold for "Fake Inclusion" removal (µm)

%% ================= MAIN PIPELINE =================

% File Selection: Matches MAPClean output format (*_clean.ctf)
fileList = dir(fullfile(dataDir, '*_clean.ctf'));
if isempty(fileList)
    % Fallback: Check for Anorthite specific files if generic clean files aren't found
    fileList = dir(fullfile(dataDir, '*_Anorthite.ctf'));
end

if isempty(fileList)
    error('GRaMC Error: No valid .ctf files found in %s. Run MAPClean first.', dataDir);
end

for fi = 1:numel(fileList)
    % --- File Setup ---
    filePath = fullfile(dataDir, fileList(fi).name);
    [~, sampleName, ~] = fileparts(fileList(fi).name);
    
    % Clean naming convention (remove common suffixes for cleaner plots)
    sampleName = erase(sampleName, '_clean'); 
    sampleName = erase(sampleName, '_Anorthite');
    
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
        
        % [Image of Twin Boundary Identification]
        % Define Twin Laws (Configured for Anorthite)
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
        S =
