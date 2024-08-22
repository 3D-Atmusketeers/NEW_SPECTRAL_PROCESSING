# Spectral-Processing

## User Guide

This guide outlines how to use the code for creating emission spectra. The process has been simplified and optimized by Isaac and Alex.

## Before Running

Before initiating the processing, ensure that the necessary data files are in place. These files are too large for GitHub storage but have been made available on Turbo.


- Download the chem data, the opacity data (including the CIA DATA), and the aerosol data. Only download what you need, these folders are quite large:
  - https://umd.app.box.com/s/klmnqdvri73wxxyylxk279bx5s6e7kxz
  - https://umd.app.box.com/s/k1g8ge5vtpji9ngn3mgiwkvyi2btzqn6
  - https://umd.app.box.com/s/32q22kh73zugvkh264msi1enf6dhr5kz
- The chemistry files should go in the DATA folder, and in a folder called chemistry_grid
- The aerosol data should go in the DATA folder, and in a folder called Aerosol_Data
- The opacity data (including the CIA opacities) should go in a folder withing DATA, called Low-Res or High-Res. Then make a folder called Full-Set and put all the opac data there.
- There should also be a Wavelenghts.txt file that you need. Ping isaacmalsky@gmail.com if you don't have it. Sorry, this is a work in progress.
- Most of the time the opacity data is too big. This will cause memory and time issues. Therefore, an extra step is needed. Run the file called crop_opacity_file.ipynb, and choose exactly what species and wavelengths you'll want. This has the benefit of dramatically reducing the size of the data files.
- Place all the General Circulation Models (GCMs) you wish to process in the `GCM-OUTPUT` directory.

## How to Run

The spectral processing suite is initiated with `run_entire_suite.py`, which can be started by executing `sbatch Run_all_sbatch` in the terminal. It usually needs at least 24 hours for runtime depending on the number of phases. Its mostly parallel. Make sure that the stuff you want to run is set correctly and the outputs are correct from run_spectra.py.

The code can also be run on a single core with run_spectra.py
This is easier than running on greatlakes, so do this unless you need a lot of models run.

The suite comprises three main steps:

1. **Altitude Regridding**: Takes about an hour per job, puts models in `/PLANET_MODELS/`.
2. **Init Files Creation**: Generates initialization files for each phase in `/Spectra/DATA/`, taking roughly 5 minutes.
3. **Final Calculations**: Outputs are saved in `/Spectra/OUT/`, with two files per phase. This step may take up to 6 hours per model

These steps are executed serially. The suite checks if altitude regridding and init files are present before proceeding, ensuring no overwrites in `/PLANET_MODELS/` or `/Spectra/DATA/`. To regenerate these files, you must manually remove or relocate the existing ones.

The suite leverages `Spectra/run_spectra.py` to perform regridding, interpolation, and post-processing. Important parameters include `NLAT` (the final number of model layers, recommend 250 or 500) and `opacity_files` (the set of opacity files in use). The `fort.7` file provides all necessary model information.

You can use these lines to only run specific wavelengths, which is helpful to double check if the code is working correctly. These values are in meters.

WAVELENGTH_START_APPROX = 0.54675e-6
WAVELENGTH_END_APPROX = 0.5468e-6
full_wavelength_range = True

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

- `Low-Res`: Medium Resolution Opacity files for general-purpose simulations.
- `High-Res`: High-resolution files for detailed analysis.

### CIA Files

The CIA (Collisionally Induced Absorption) files are crucial for providing the opacity data related to molecular interactions that are not captured by line transitions alone. These files account for the absorption caused by collisions between atmospheric molecules, contributing significantly to the overall opacity in dense, cooler parts of an atmosphere.

### Chem Files

Chem files detail the atmospheric abundances of various molecules and elements. Chem files are shared between opacity sets. These files are located in a common directory accessible to all opacity sets. There are a number of chem files available based on metallicity. If your naming conventions are correct, then it should automatically select the closet one.

### Rayleigh Scattering

In addition to the opacities provided by the EOS, CIA, and line data files, Rayleigh scattering is another source of opacity considered by the suite.

When updating or adding new data sets, please ensure this README is kept current with the latest file organization and opacity set descriptions.


## Notes

Running the suite generates numerous files in `/Spectra/`, typical for parallel phase processing. To clean up, execute `Clean_suite.py` to remove unnecessary files. Be aware that numerous Slurm files will also be created.

Currently, the suite cannot process multiple planetary systems simultaneously due to varying stellar parameters required for each run. `run_spectra` allows for single-core planet processing, whereas `set_up_spectra_folders.py` enables concurrent processing of multiple planets, though this method involves significant data duplication.

### Required Files

- `fort.26`, `fort.50`, `fort.7`, `fort.62`

**Important**: `set_up_spectra_folders` automates folder and file management but requires post-run reorganization to ensure proper file placement.
