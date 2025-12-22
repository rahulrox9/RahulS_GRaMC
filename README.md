# GRaMC – Grain Reconstruction and Multi-stage Cleaning

**Note:** GRaMC is designed to be used in tandem with **MAPClean** (Microstructurally Adaptive Pixel-Level Cleaning), the primary software associated with Subbaraman et al. (in prep). While **MAPClean** handles pixel-level noise removal and data restoration, **GRaMC** takes the clean output and handles **grain-level logic** to produce robust grain sets for analysis.

It moves beyond standard MTEX grain reconstruction by incorporating microstructural logic specifically for complex geological materials (e.g., Plagioclase). It handles twin boundary merging, fake inclusion removal, and matrix isolation to produce grains suitable for high-level shape and texture analysis.

## Key Features
* **Multi-Stage Reconstruction:** Checkpointed workflow allowing users to inspect grains at every stage of the reconstruction.
* **Twin Boundary Merging:** Automatically detects and merges grains separated by specific twin laws. Currently pre-configured for **Anorthite** laws (Albite, Pericline, Carlsbad, Manebach, Baveno).
* **Fake Inclusion Removal:** Identifies small or zero-GOS "grains" trapped inside larger host grains and merges them into the host, cleaning up the microstructure without losing data.
* **Background Isolation:** Filters grain populations by size to isolate large phenocrysts from the groundmass/matrix.
* **Advanced Property Calculation:** Automatically recalculates **Grain Orientation Spread (GOS)** and **Mean Orientations** during merging events using area-weighted logic.

## Requirements
* **MATLAB** (Tested on 2024b)
* **MTEX Toolbox** (Tested on v 6.0.0)
* **Image Processing Toolbox**

## Installation
1. Clone or download the repository (often located within the MAPClean suite).
2. Ensure the **MAPClean** outputs (cleaned `.ctf` files) are located in the `DataFiles` folder.
3. Add the script to your MATLAB path.

## Usage
1. Open `GRaMC.m`.
2. Set the **Stage Control Flags** to determine which steps to run:
```matlab
run_initial   = false;  % Initial reconstruction
run_twinMerge = true;   % Merge twins
run_fakeInc   = true;   % Remove fake inclusions
run_bgRemoval = true;   % Isolate large grains
run_Final     = true;   % Save final output
```
3. Run the script.
4. Final output files (`*_finalGrains.mat`) are saved in the `checkpoints` directory, ready for downstream analysis (e.g. using GRaFT for fabric and texture analysis).

## Workflow Details
1. **Initial Reconstruction:**  Grains are calculated from pixel data using a misorientation threshold (default 10°). Small noise grains are immediately removed.
2. **Twin Merging:** Child grains sharing a specific twin boundary (e.g., Albite law) are merged into a single parent grain to reconstruct the original growth morphology.
3. **Inclusion Cleaning:** Small "islands" inside larger grains (often artifacts of pixel cleaning) are absorbed into the host grain.
4. **Matrix Removal:** A size filter isolates the primary grain population (e.g., macrocrysts) for shape analysis.

## Parameters

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `grainThreshold` | 10° | Misorientation angle threshold for defining grain boundaries. |
| `minPixelsGrains` | 5 | Minimum number of pixels required to constitute a valid grain. |
| `fakeIncSizeThres`| 100 µm | Maximum size for an "inclusion" to be merged into its host. |
| `bgSizeThres` | 300 µm | Size threshold for Background Removal. Grains smaller than this are discarded in the final step (useful for separating phenocrysts from matrix). |
| `exportRes` | 300 dpi | Resolution for exported PNG maps. |

## Directory Structure
This tool expects to run downstream of MAPClean.
```text
ProjectRoot/
├── MAPClean.m          # (Step 1) Pixel Cleaning
├── GRaMC.m             # (Step 2) Grain Reconstruction (This Tool)
├── DataFiles/          # Input .ctf files (Exported from MAPClean)
├── checkpoints/        # Intermediate .mat saves (initial, twin, etc.)
├── exports/
│   └── GrainClean/
│       └── [SampleID]/ # Diagnostic PNG maps (Boundaries, GOS, IPF)

└── README.md
```
## License
This code is licensed under **GPL version 3**.
