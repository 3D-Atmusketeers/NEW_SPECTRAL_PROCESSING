#!/usr/bin/env python
import numpy as np
from sys import exit
import os
import run_grid
import altitude_regridding
import add_clouds
import convert_fort_files
import re
import subprocess
import grab_input_data
from setup_opac_versions import replace_files, modify_input_h, insert_opacity_definitions, modify_totalopac
from get_wavelength_grid import find_closest_wavelength_indices
from chemistry_file_selector import find_closest_chemistry_file
#import Clean_suite

# Phases in degrees, inclination in radians (sorry)
# An inclination of 0 corresponds to edge on
phases = [0.0]
inclinations = [0.0]
system_obliquity = 0

# I recommend leaving these as is
# The NLAT and NLON can be changed, but these values work well
NTAU = 250

# Cut of the bottom of the atmosphere if needed
# DONT MESS WITH THIS, ISAAC HASN'T FULLY CODED IT!!!!!
max_pressure_bar = 100

# Please don't touch these
NLAT = 48
NLON = 96

# 0 is off
# 1 is everything
# 2 is Wind only
# 3 is rotation only
dopplers = [0]

# If you only need to change the phase you can use this knob
# It skips a lot of steps for the regridding
ONLY_PHASE = True

# If you only have the fort files use this
# Please don't only have the fort files
# It requires that the fort files are named according to something particular
# In the correct directory
USE_FORT_FILES = True

# This is a polynomial fit of the GCM TP profiles in order to smooth it for the post processing
# If false, just use the GCM TP profiles
smoothing = True

# These are the planet files that you need to run the code
# They should be pretty big files, and don't include the .txt with the names here
planet_names = ["HD189-DOGRAY"]

# The options are lowres and hires
# Isaac Malsky is still working on highres
opacity_set_id = 'Low-Res'

# Specify the wavelength range that you'd like to calculate
# If values aren't given, or if they're negative -1 for both
# Then it will calculate the entire grid
WAVELENGTH_START_APPROX = 0.54675e-6
WAVELENGTH_END_APPROX = 0.5468e-6
full_wavelength_range = True
LAMBDA_START, LAMBDA_END, START_WAVELENGTH, END_WAVELENGTH, NLAMBDA = find_closest_wavelength_indices(opacity_set_id,
                                                                                             full_wavelength_range,
                                                                                             WAVELENGTH_START_APPROX,
                                                                                             WAVELENGTH_END_APPROX)

# Construct the path to the directory
opacity_files_directory = os.path.join('DATA', opacity_set_id)

# Adjust the list comprehension to parse filenames
opacity_species = [file[4:-4] for file in os.listdir(opacity_files_directory)
                   if file.startswith("opac") and "CIA" not in file and file.endswith(".dat")]

# Check if H2O, CO, and CO2 are included
required_species = ["H2O"]
missing_species = [species for species in required_species if species not in opacity_species]

if missing_species:
    print("\n" + "="*60)
    print("WARNING: The following species are not included:", ', '.join(missing_species))
    print("Stopping the code to prevent incorrect calculations.")
    print("="*60 + "\n")
    # Stop the code here to prevent further execution
    quit()

##################################################
######        SET THE OPAC SPECIES        ########
##################################################
opacity_species = ['H2O']

print("\n" + "="*60)
print("WARNING: Using a limited subset of available species!")
print("-" * 60)
print(f"Selected Species: {opacity_species}")
print("="*60 + "\n")

# Set the wavelength to evaluate the clouds at for plotting!
# This could be put in a better place I think
cloud_file = 'DATA/Aerosol_Data/wavelength_array_for_cloud_scattering_data_in_microns.txt'
if os.path.exists(cloud_file):
    wavelength_grid = np.loadtxt(cloud_file)
else:
    print("Warning: You're missing the cloud files that should be in DATA/Aerosol_Data")
    exit()

if opacity_set_id == 'Low-Res':
    cloud_wavelength = 0.500
    wav_loc = np.absolute(wavelength_grid-cloud_wavelength).argmin()
else:
    cloud_wavelength = 2.00
    wav_loc = np.absolute(wavelength_grid-cloud_wavelength).argmin()

for q in range(len(planet_names)):
    print()
    planet_name = planet_names[q]

    # Assuming planet_names and q are defined somewhere above this
    planet_name_original = planet_names[q]

    # Get the current working directory
    current_directory = os.getcwd()

    runname     = planet_name + '/Planet_Run'
    path        = '../GCM-OUTPUT/'

    aerosol_layers = int(grab_input_data.get_input_data(path, runname,"fort.7", "AERLAYERS"))
    grav           = grab_input_data.get_input_data(path, runname, "fort.7","GA")
    gasconst       = grab_input_data.get_input_data(path, runname, "fort.7","GASCON")
    R_PLANET       = grab_input_data.get_input_data(path, runname, "fort.7","RADEA")
    P_ROT          = (grab_input_data.get_input_data(path, runname, "fort.7","WW") / (2.0*np.pi)*(24*3600)) ** -1.0
    oom            = grab_input_data.get_input_data(path, runname, "fort.7","OOM_IN")
    MTLX           = grab_input_data.get_input_data(path, runname, "fort.7","MTLX")
    MET_X_SOLAR    = 10.0 ** grab_input_data.get_input_data(path, runname, "fort.7","METALLICITY")
    HAZES = grab_input_data.get_input_data(path, runname, "fort.7", "HAZES")[0] == 'T'
    MOLEF          = grab_input_data.get_input_data(path, runname, "fort.7","MOLEF")
    HAZE_TYPE = next((s for s in ['soot', 'soot_2xpi0', 'sulfur', 'tholin'] if s in runname.lower()), 'None')


    GAS_CONSTANT_R = 8.314462618
    MEAN_MOLECULAR_WEIGHT = np.round((GAS_CONSTANT_R/gasconst) * 1000, 4)

    # This is the path to the chemistry file
    chemistry_file_path = find_closest_chemistry_file(MET_X_SOLAR)

    print("\n" + "="*40)
    print("======== RUNNING SIMULATION ========")
    print(f"Planet Name: {planet_name}")
    print(f"NLAMBDA: {NLAMBDA}")
    print(f"STARTING WAVELENGTH: {START_WAVELENGTH}")
    print(f"ENDING WAVELENGTH: {END_WAVELENGTH}")
    print("="*40 + "\n")

    # Check if the substring 'soot_2xpi0' is in HAZE_TYPE
    if 'soot_2xpi0' in HAZE_TYPE:
        # Replace 'soot_2xpi0' with 'soot-2xpi0'
        HAZE_TYPE = HAZE_TYPE.replace('soot_2xpi0', 'soot-2xpi0')

    # Whether  there are clouds
    # 0 is no clouds, 1 is clouds
    # This is also important for filling in the correct number of 0s for the input files
    if any(i > 1e-20 for i in MOLEF):
        CLOUDS = 1
    else:
        CLOUDS = 0

    # These values are used mostly for the fort files
    # A couple of them are also used for the spectra
    # So make sure all the constants are set correctly
    # Eventually these should be pulled from the fort.7 file
    # Alas, I am lazy, so you have to do it by hand
    if os.path.isfile(path+runname+'/fort.26'):
        with open(path + runname + '/fort.26') as f:
            lines_fort26 = f.readlines()
            INITIAL_NTAU = int(float(re.findall(r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?", lines_fort26[0])[2]))
    elif os.path.isfile(path+runname+'/fort.2600'):
        with open(path + runname + '/fort.2600') as f:
            lines_fort26 = f.readlines()
            INITIAL_NTAU = int(float(re.findall(r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?", lines_fort26[0])[2]))
    else:
        print()
        print("!!!!!!!!!! You need the fort.26 or fort.2600 file !!!!!!!!!!!")
        print()
        exit(0)

    # get all the stellar parameters needed from the excel doc of stellar params
    # that doc should be in the spectra folder, make sure to update it every once in a while
    star_name    = grab_input_data.read_planet_and_star_params(planet_name, "Name")
    ORB_SEP      = float(grab_input_data.read_planet_and_star_params(planet_name, "a (au)")) * 1.496e11
    STELLAR_TEMP = float(grab_input_data.read_planet_and_star_params(planet_name, "T* (K)"))
    R_STAR       = float(grab_input_data.read_planet_and_star_params(planet_name, "R* (R_sun)")) * 695700000

    print("\n" + "="*60)
    print("Planet Characteristics")
    print("\tCloud Types in Order with Corresponding Amounts:")
    cloud_types = "KCl, ZnS, Na2S, MnS, Cr, SiO2, Mg2SiO4, VO, Ni, Fe, Ca2SiO4, CaTiO3, Al2O3"
    print(f"\tCloud Types: {cloud_types}")
    print(f"\tMOLEF: {MOLEF}")
    print(f"\tGravity: {grav} m/s^2")
    print(f"\tNumber of Orders of Magnitude (OOM): {oom}")
    print(f"\tGas Constant: {gasconst} J/(K*mol)")
    print(f"\tPlanet Radius: {R_PLANET} km")
    print(f"\tOrbital Period: {P_ROT} days")
    print(f"\tMetallicity (MTLX): {MTLX}")
    print(f"\tHazes Present: {'Yes' if HAZES else 'No'}")
    print(f"\t{aerosol_layers} Cloud Layers")
    print(f"\tHaze Type: {HAZE_TYPE}")
    print(f"\tGCM Layers: {INITIAL_NTAU}")
    print(f"\tMean Molecular Weight: {MEAN_MOLECULAR_WEIGHT}")
    print(f"\tCloud wavelength (microns): {cloud_wavelength}")

    print("\nOpacity and Chemistry information")
    print(f"\tMetallicity relative to solar (METALLICITY): {MET_X_SOLAR}")
    print(f"\tChemistry File Path: {chemistry_file_path}")
    print(f"\tUsing Opacity Set: {opacity_set_id}")
    print("="*60 + "\n")

    print("\n" + "="*60)
    print("Star Characteristics")
    print(f"\tYou matched '{star_name}' on planet '{planet_name}'")
    print(f"\tOrbital Separation: {ORB_SEP / 1.496e11:.2f} AU")
    print(f"\tStar Temp: {STELLAR_TEMP} K")
    print(f"\tStar Radius: {R_STAR / 695700000:.2f} Solar Radii")
    print("="*60 + "\n\n")


    # Are these used?
    surfp = 100 #surface pressure, in bars
    tgr   = 1500 #temperature at 100 bars

    # This is specifically for the regridding
    # You can change the grid density
    # I wouldn't mess with this though
    grid_lat_min = -87.16
    grid_lat_max = 87.16
    grid_lon_min = 0.0
    grid_lon_max = 356.25

    def get_run_lists(phases, inclinations):
        for phase in phases:
            for inc in inclinations:
                phase = str(phase)
                inc = str(inc)
                input_paths.append('DATA/init_' + planet_name + '_phase_{}_inc_{}.txt'.format(phase, inc))
                inclination_strs.append(inc)
                phase_strs.append(phase)

        return input_paths, inclination_strs, phase_strs

    def run_exo(input_paths, inclination_strs, phase_strs, doppler_val):
        """
        This function orchestrates running Eliza's code for each provided input path,
        modifying 'input.h' appropriately for each run based on dynamic parameters.
        """
        # Prepare the base files (e.g., copy templates without 'template_' prefix)
        replace_files(opacity_set_id)

        # Generate output paths based on input paths and doppler value
        output_paths = ['OUT/Spec_' + str(doppler_val) + '_' + path[10:-4] for path in input_paths]

        for i, input_path in enumerate(input_paths):
            # Construct modifications dictionary for this run
            modifications = {
                "<<output_file>>": "\"" + output_paths[i] + "\"",
                "<<input_file>>": "\"" + input_path + "\"",
                "<<doppler>>": str(doppler_val),
                "<<inclination>>": inclination_strs[i],
                "<<phase>>": phase_strs[i],
                "<<CHEMISTRY_FILE>>": "\"" + chemistry_file_path + "\"",
                "<<CLOUDS>>":str(CLOUDS),
                "<<NTAU>>":str(NTAU),
                "<<NLAT>>":str(NLAT),
                "<<NLON>>":str(NLON),
                "<<NLAMBDA>>":str(NLAMBDA),
                "<<LAMBDA_START>>":str(LAMBDA_START),
                "<<LAMBDA_END>>":str(LAMBDA_END),
                "<<W0_VAL>>":str(W0_VAL),
                "<<G0_VAL>>":str(G0_VAL),
                "<<GRAVITY_SI>>":str(grav),
                "<<R_PLANET>>":str(R_PLANET),
                "<<ORB_SEP>>":str(ORB_SEP),
                "<<STELLAR_TEMP>>":str(STELLAR_TEMP),
                "<<R_STAR>>":str(R_STAR),
                "<<P_ROT>>":str(P_ROT),
                "<<MEAN_MOLECULAR_WEIGHT>>":str(MEAN_MOLECULAR_WEIGHT),
                "<<HAZE_TYPE>>":"\"" + HAZE_TYPE +"\"",
                "<<HAZES>>": str(1) if HAZES else str(0),
                }

            # Call the function to modify 'input.h' for this run
            modify_input_h(modifications, opacity_set_id)

            # Update the input.h file for the species desired
            insert_opacity_definitions('input.h', 'DATA/' + opacity_set_id, opacity_species)

            # Modify totalopac.c
            modify_totalopac(opacity_species)

            try:
                # Run Eliza's code
                subprocess.run('make rt_emission_aerosols.exe', shell=True, check=True)
                file_name = f"rt_emission_aerosols_{planet_name}_phase_{phase_strs[i]}.exe"
                os.rename('rt_emission_aerosols.exe', file_name)

                permission_command = f'chmod 755 {file_name}'
                subprocess.run(permission_command, shell=True, check=True)

                execution_command = f"./{file_name}"
                subprocess.run(execution_command, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                # Throw an error
                print(f"An error occurred while executing a subprocess: {e}")
        return None


    input_paths = []
    output_paths = []
    inclination_strs = []
    phase_strs = []

    STEP_ONE = False
    STEP_TWO = False
    STEP_THREE = True

    if STEP_ONE:
        # Convert the fort files to the correct format

        if USE_FORT_FILES == True:
            convert_fort_files.convert_to_correct_format(path, runname, planet_name, INITIAL_NTAU, surfp, oom, tgr, grav, gasconst)
            print ("Converted the fort files to the new format")
        else:
            pass

        add_clouds.add_clouds_to_gcm_output(path, runname, planet_name,
                                            grav, MTLX, CLOUDS, MOLEF,
                                            aerosol_layers, INITIAL_NTAU,
                                            gasconst, HAZE_TYPE, HAZES, wav_loc, MET_X_SOLAR)


        # Regrid the file to constant altitude and the correct number of layers
        altitude_regridding.regrid_gcm_to_constant_alt(path, CLOUDS, planet_name, NLAT, NLON, INITIAL_NTAU, NLON, NTAU, HAZES, max_pressure_bar, smoothing)

        print ("Regridded the planet to constant altitude")
    if STEP_TWO:
        # If you already have the Final planet file creates you can commend out run_grid and double planet file
        run_grid.run_all_grid(planet_name, phases, inclinations, system_obliquity, NTAU, NLAT, NLON, grid_lat_min, grid_lat_max, grid_lon_min, grid_lon_max, ONLY_PHASE)
    if STEP_THREE:
        # Get all the files that you want to run
        input_paths, inclination_strs, phase_strs = get_run_lists(phases, inclinations)

        # If you want to manually set these values you can leave them here
        # Normally they will not affect it, unless you manually set them in two_stream.h
        W0_VALS = [0.0]
        G0_VALS = [0.0]

        for G0_VAL in G0_VALS:
            for W0_VAL in W0_VALS:
                for doppler_val in dopplers:
                    run_exo(input_paths, inclination_strs, phase_strs, doppler_val)


    #print('Finished running', planet_name)


def cleanup_lock_and_exe_files():
    """
    Cleans up any files with 'lock' in their name or ending with '.exe' in the current directory.
    """
    try:
        files_to_remove = [f for f in os.listdir('.') if 'lock' in f or f.endswith('.exe')]
        for file in files_to_remove:
            os.remove(file)
            #print(f"Removed file: {file}")
    except Exception as e:
        print(f"An error occurred while cleaning up files: {e}")

# Call the function to clean up lock files
#cleanup_lock_and_exe_files()

#uncomment this out if you would like the files to automatically delete the bonus files that are created
#Clean_suite.automaticclean(__file__)
