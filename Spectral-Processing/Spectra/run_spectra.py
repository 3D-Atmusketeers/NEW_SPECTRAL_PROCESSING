#!/usr/bin/env python

from shutil import copyfile
from sys import exit
import os
import sys
import numpy as np
import run_grid
import altitude_regridding
import add_clouds
import convert_fort_files
import re
import shutil
import grab_input_data
import setup_opac_versions

# Phases in degrees, inclination in radians (sorry)
# An inclination of 0 corresponds to edge on
phases = [0.0]
inclinations = [0.0]
system_obliquity = 0

# I recommend leaving these as is
# The NLAT and NLON can be changed, but these values work well
NTAU = 250

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

# There are the different sets of opacity and EOS files
# There are somethings that need to be changed in the template inputs file to make this happen
# If you change the underlying data files these might need to be changed
opacity_files = 'SET_1'

# These are the planet files that you neesd to run the code
# They should be pretty big files, and don't include the .txt with the names here
planet_names = ["GJ1214b-CLEAR-100X"]

# Set the wavelength to evaluate the clouds at for plotting!
# This could be put in a better place I think
wavelength_grid = np.loadtxt('SCATTERING_DATA/wavelength_array_for_cloud_scattering_data_in_microns.txt')
if opacity_files == 'SET_1':
    cloud_wavelength = 5.0
    wav_loc = np.absolute(wavelength_grid-cloud_wavelength).argmin()
elif opacity_files == 'SET_2':
    cloud_wavelength = 2.3
    wav_loc = np.absolute(wavelength_grid-cloud_wavelength).argmin()
elif opacity_files == 'SET_3':
    cloud_wavelength = 2.3
    wav_loc = np.absolute(wavelength_grid-cloud_wavelength).argmin()
else:
    print("YOU NEED TO SET WHICH OPACITY SET YOU'RE USING")
    exit(0)

for q in range(len(planet_names)):
    planet_name = planet_names[q]


    runname     = planet_name + '/Planet_Run'
    path        = '../GCM-OUTPUT/'

    MOLEF          = grab_input_data.get_input_data(path, runname, "fort.7", "MOLEF")
    aerosol_layers = int(grab_input_data.get_input_data(path, runname,"fort.7", "AERLAYERS"))
    grav           = grab_input_data.get_input_data(path, runname, "fort.7","GA")
    gasconst       = grab_input_data.get_input_data(path, runname, "fort.7","GASCON")
    R_PLANET       = grab_input_data.get_input_data(path, runname, "fort.7","RADEA")
    P_ROT          = (grab_input_data.get_input_data(path, runname, "fort.7","WW") / (2.0 * np.pi) * (24 * 3600)) ** -1.0
    oom            = grab_input_data.get_input_data(path, runname, "fort.7","OOM_IN")
    MTLX           = grab_input_data.get_input_data(path, runname, "fort.7","MTLX")
    HAZES          = grab_input_data.get_input_data(path, runname, "fort.7","HAZES")

    # Necessary for choosing the chem table!
    MET_X_SOLAR    = 10.0 ** grab_input_data.get_input_data(path, runname, "fort.7","METALLICITY")

    if (opacity_files == "SET_1"):
        if (0.9 * MET_X_SOLAR <= 1.0 <= 1.1 * MET_X_SOLAR):
            chemistry_file_path = "DATA/SET_1/ordered_1x_solar_metallicity_chem.dat"
        elif (0.9 * MET_X_SOLAR <= 100.0 <= 1.1 * MET_X_SOLAR):
            chemistry_file_path = "DATA/SET_1/ordered_100x_solar_metallicity_chem.dat"
        elif (0.9 * MET_X_SOLAR <= 300.0 <= 1.1 * MET_X_SOLAR):
            chemistry_file_path = "DATA/SET_1/ordered_300x_solar_metallicity_chem.dat"
        elif (0.9 * MET_X_SOLAR <= 3000.0 <= 1.1 * MET_X_SOLAR):
            chemistry_file_path = "DATA/SET_1/ordered_3000x_solar_metallicity_chem.dat"
        else:
            print("Error in choosing which metallicy the chemistry file should be")

    elif (opacity_files == "SET_2"):
        chemistry_file_path = "DATA/SET_2/eos_solar_doppler.dat"
    elif (opacity_files == "SET_3"):
        chemistry_file_path = "DATA/SET_3/eos_solar_doppler.dat"
    else:
        print("Error in choosing the chemistry file!")

    GAS_CONSTANT_R = 8.314462618
    GASCON = grab_input_data.get_input_data(path, runname, "fort.7","GASCON")
    MEAN_MOLECULAR_WEIGHT = np.round((GAS_CONSTANT_R/GASCON) * 1000, 4)

    print ()
    print ()
    print ("*************************************")
    print ("*************************************")
    print ("RUNNING: " + planet_name)

    # THIS IS A PROBLEM BECAUSE SOOT IS NOT DEFINED IN fort.7
    if ("SOOT".lower() in planet_names[q].lower()):
        HAZE_TYPE = 'soot'
    elif ("THOLIN".lower() in planet_names[q].lower()):
        HAZE_TYPE = 'tholin'
    else:
        HAZE_TYPE = 'None'

    if (HAZE_TYPE == 'soot' or HAZE_TYPE == 'soot-2xpi0' or HAZE_TYPE == 'tholin' or HAZE_TYPE == 'sulfur'):
        HAZES = True
    else:
        HAZES = False

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
    else:
        print ("You need the fort.26 file")

    # get all the stellar parameters needed from the excel doc of stellar params
    # that doc should be in the spectra folder, make sure to update it every once in a while
    ORB_SEP      = grab_input_data.read_planet_and_star_params(planet_name, "a (au)") * 1.496e11
    STELLAR_TEMP = grab_input_data.read_planet_and_star_params(planet_name, "T* (K)")
    R_STAR       = grab_input_data.read_planet_and_star_params(planet_name, "R* (R_sun)") * 695700000

    print("Planet characteristics")
    print("These are the cloud types in order, and the corresponding amounts")
    print ("KCl, ZnS, Na2S, MnS, Cr, SiO2, Mg2SiO4, VO, Ni, Fe, Ca2SiO4, CaTiO3, Al2O3")
    print("MOLEF", MOLEF)
    print("Gravity =  ",     grav)
    print("Number of OOM = ", oom)
    print("Gas Constant = ", gasconst)
    print("Planet Radius = ", R_PLANET)
    print("Orbital Period = ", P_ROT)
    print("MTLX = ", MTLX)
    print("There are hazes? ", HAZES)
    print(aerosol_layers, "aerosol layers")
    print("Haze type:", HAZE_TYPE)
    print("GCM Layers = ", INITIAL_NTAU)
    print("Mean Molecular Weight", MEAN_MOLECULAR_WEIGHT)
    print("")
    print("Be careful to make sure that your chemistry file is correct!")
    print("METALLICITY = ", MET_X_SOLAR, chemistry_file_path)

    print("")
    print("Star characteristics")
    print("Orbital Separation = ",ORB_SEP)
    print("Star Temp = ",STELLAR_TEMP)
    print("Star Radius = ",R_STAR)


    # Are these used?
    surfp = 100 #surface pressure, in bars
    tgr   = 1500 #temperature at 100 bars

    print ("*************************************")
    print ("*************************************")

    print ("")
    print ("")

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
        This runs Eliza's code
        """
        inputs_file = 'input.h'
        output_paths = []

        # The output paths should be similar to the input paths
        # Minus the .dat file extension and saved to OUT/
        for file_path in input_paths:
            output_paths.append('OUT/Spec_' + str(doppler_val) + '_' + file_path[10:-4])

        # Each Run needs to have a specific input.h file
        # With the correct input and output paths
        for i in range(len(input_paths)):
            output_temp = output_paths[i]
            input_temp  = input_paths[i]
            
            # This will copy all the files that need to be changed over with different opac versions
            setup_opac_versions.replace_files(opacity_files, MET_X_SOLAR)

            # Copy the template for inputs
            try:
                copyfile('template_inputs.h', inputs_file)
            except IOError as e:
                print("Unable to copy file. %s" % e)
                exit(1)
            except:
                print("Unexpected error:", sys.exc_info())
                exit(1)

            # Read in the file
            with open(inputs_file, 'r') as file :
                filedata = file.read()

            # Replace the input and output paths
            filedata = filedata.replace("<<output_file>>", "\"" + output_temp + str(W0_VAL) + str(G0_VAL) + "\"")
            filedata = filedata.replace("<<input_file>>", "\"" + input_temp + "\"")
            filedata = filedata.replace("<<doppler>>", str(doppler_val))
            filedata = filedata.replace("<<inclination>>", inclination_strs[i])
            filedata = filedata.replace("<<phase>>", phase_strs[i])
            filedata = filedata.replace("<<CHEMISTRY_FILE>>", "\"" + chemistry_file_path + "\"")

            filedata = filedata.replace("<<CLOUDS>>", str(CLOUDS))
            filedata = filedata.replace("<<NTAU>>", str(NTAU))
            filedata = filedata.replace("<<NLAT>>", str(NLAT))
            filedata = filedata.replace("<<NLON>>", str(NLON))

            filedata = filedata.replace("<<W0_VAL>>", str(W0_VAL))
            filedata = filedata.replace("<<G0_VAL>>", str(G0_VAL))

            filedata = filedata.replace("<<GRAVITY_SI>>", str(grav))

            filedata = filedata.replace("<<R_PLANET>>",     str(R_PLANET))
            filedata = filedata.replace("<<ORB_SEP>>",      str(ORB_SEP))
            filedata = filedata.replace("<<STELLAR_TEMP>>", str(STELLAR_TEMP))
            filedata = filedata.replace("<<R_STAR>>",       str(R_STAR))
            filedata = filedata.replace("<<P_ROT>>",        str(P_ROT))

            filedata = filedata.replace("<<MEAN_MOLECULAR_WEIGHT>>", str(MEAN_MOLECULAR_WEIGHT))


            filedata = filedata.replace("<<HAZE_TYPE>>", "\"" + HAZE_TYPE +"\"")
            if (HAZES == True):
                filedata = filedata.replace("<<HAZES>>", str(1))
            else:
                filedata = filedata.replace("<<HAZES>>", str(0))

            # Write the file out again
            with open(inputs_file, 'w') as file:
                file.write(filedata)
            
            # Run Eliza's code
            os.system('make clean')
            os.system('make rt_emission_aerosols.exe')
            os.system('./rt_emission_aerosols.exe')

        return None


    input_paths = []
    output_paths = []
    inclination_strs = []
    phase_strs = []


    # Convert the fort files to the correct format
    if USE_FORT_FILES == True:
        convert_fort_files.convert_to_correct_format(path, runname, planet_name, INITIAL_NTAU, surfp, oom, tgr, grav, gasconst)
        print ("Converted the fort files to the new format")
    else:
        pass

    add_clouds.add_clouds_to_gcm_output(path, runname, planet_name,
                                        grav, MTLX, CLOUDS, MOLEF,
                                        aerosol_layers, INITIAL_NTAU,
                                        gasconst, HAZE_TYPE, HAZES, wav_loc)

    # Regrid the file to constant altitude and the correct number of layers
    altitude_regridding.regrid_gcm_to_constant_alt(path, CLOUDS, planet_name, NLAT, NLON, INITIAL_NTAU, NLON, NTAU, HAZES)
    print ("Regridded the planet to constant altitude")

    # If you already have the Final planet file creates you can commend out run_grid and double planet file
    run_grid.run_all_grid(planet_name, phases, inclinations, system_obliquity, NTAU, NLAT, NLON, grid_lat_min, grid_lat_max, grid_lon_min, grid_lon_max, ONLY_PHASE)

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

print("Moving the files out of the clean directory")
for filename in os.listdir('DATA'):
    if re.match(r'init_*', filename):
        shutil.move(os.path.join('DATA', filename), os.path.join('../PLANET_MODELS', filename))

for filename in os.listdir('OUT'):
    if re.match(r'Spec_*', filename):
        shutil.move(os.path.join('OUT', filename), os.path.join('../FINISHED_SPECTRA', filename))