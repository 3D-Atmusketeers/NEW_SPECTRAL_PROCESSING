# Spectral-Processing

## User Guide

This guide outlines how to use the code for creating emission spectra. The process has been simplified and optimized by Isaac and Alex.

## Before Running

Before initiating the processing, ensure that the necessary data files are in place. These files are too large for GitHub storage but have been made available on Turbo.

- Run `download_all_data_files_from_turbo.py` to download the required files. Check the `files_to_download` list at the bottom of this script to see which files are being retrieved.
- Once downloaded, move the `SET_X` data folder to `Spectra/DATA`, and the `SCATTERING_DATA` to `Spectra/SCATTERING_DATA`.
- Place all the General Circulation Models (GCMs) you wish to process in the `GCM-OUTPUT` directory.
- Set 1 should be medium resolution, and Set 2 should be high resolution

## How to Run

The spectral processing suite is initiated with `run_entire_suite.py`, which can be started by executing `sbatch Run_all_sbatch` in the terminal. It usually needs at least 24 hours for runtime depending on the number of phases. Its mostly parallel.

The code can also be run on a single core with run_spectra.py
This is easier than running on greatlakes, so do this unless you need a lot of models run.

The suite comprises three main steps:

1. **Altitude Regridding**: Takes about an hour per job, puts models in `/PLANET_MODELS/`.
2. **Init Files Creation**: Generates initialization files for each phase in `/Spectra/DATA/`, taking roughly 5 minutes.
3. **Final Calculations**: Outputs are saved in `/Spectra/OUT/`, with two files per phase. This step may take up to 6 hours per model

These steps are executed serially. The suite checks if altitude regridding and init files are present before proceeding, ensuring no overwrites in `/PLANET_MODELS/` or `/Spectra/DATA/`. To regenerate these files, you must manually remove or relocate the existing ones.

The suite leverages `Spectra/run_spectra.py` to perform regridding, interpolation, and post-processing. Important parameters include `NLAT` (the final number of model layers, recommend 250 or 500) and `opacity_files` (the set of opacity files in use). The `fort.7` file provides all necessary model information.

For cloud-dependent wavelength properties adjustments, modifications in `run_spectra.py` are necessary.

As a user, PLEASE CHECK THAT THE OUTPUT VARIABLES ARE CORRECT!!!!

## Visualizations

Run `create_all_figures.py` to generate all coded figures. Post-processing may leave empty folders, which should be cleaned up manually.

## Folders Overview

- `Figures`: Contains generated figures.
- `Visualizations`: Stores visualization code.
- `FINISHED_SPECTRA`: Destination for completed spectral outputs.
- `GCM-OUTPUT`: Input folder for GCM outputs.
- `PLANET_MODELS`: Stores regridded GCMs and interpolation results.
- `Spectra`: Main processing directory; avoid placing raw data here.

## File Naming

The post-processing automatically extracts necessary information from `fort.7`, removing the need for specific naming conventions for GCM outputs. Ensure different planet titles to prevent data overwrites. The output printed to the terminal provides a summary of processed data.

## Opacity Versions and Additional Files

The suite supports various resolution and temperature regimes, necessitating appropriate EOS, opacity files, and additional data files for comprehensive atmospheric modeling. While swapping out these files for different simulations is straightforward, extensive modifications to the codebase to accommodate new versions are not recommended. The relevant files are organized within the `OPAC_CODE_VERSIONS` and `DATA` directories, and further distinctions are made between different types of files:

- `SET_1`: Medium Resolution Opacity files for general-purpose simulations.
- `SET_3`: High-resolution files for detailed analysis.

### CIA Files

The CIA (Collisionally Induced Absorption) files are crucial for providing the opacity data related to molecular interactions that are not captured by line transitions alone. These files account for the absorption caused by collisions between atmospheric molecules, contributing significantly to the overall opacity in dense, cooler parts of an atmosphere. CIA files are specific to each set of opacity files and reside within the corresponding `SET_X` directories under `OPAC_CODE_VERSIONS`.

### Chem Files

Chem files detail the atmospheric abundances of various molecules and elements. Unlike the CIA files, chem files are shared between opacity sets. These files are located in a common directory accessible to all opacity sets.

### Rayleigh Scattering

In addition to the opacities provided by the EOS, CIA, and line data files, Rayleigh scattering is another source of opacity considered by the suite.

When updating or adding new data sets, please ensure this README is kept current with the latest file organization and opacity set descriptions.


## Notes

Running the suite generates numerous files in `/Spectra/`, typical for parallel phase processing. To clean up, execute `Clean_suite.py` to remove unnecessary files. Be aware that numerous Slurm files will also be created.

Currently, the suite cannot process multiple planetary systems simultaneously due to varying stellar parameters required for each run. `run_spectra` allows for single-core planet processing, whereas `set_up_spectra_folders.py` enables concurrent processing of multiple planets, though this method involves significant data duplication.

### Required Files

- `fort.26`, `fort.50`, `fort.7`, `fort.62`

**Important**: `set_up_spectra_folders` automates folder and file management but requires post-run reorganization to ensure proper file placement.
