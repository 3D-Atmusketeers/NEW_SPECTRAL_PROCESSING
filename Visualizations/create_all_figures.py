#!/usr/bin/python
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import re

# These are the plotting codes that Isaac wrote
# If they're broken, venmo isaac-malsky
# I wont' fix them, but I'll appreciate the money
import broadband_phase_curves
import cloud_coverage_isobars
import pressure_temperature_condensation_curves
import spectra

import grab_input_data

import aerosol_maps
import aerosol_profiles
import wind_maps
#import emission_maps

# Make it pretty!
plt.style.use('science.mplstyle')


# Figure out what planets!
planet_names = [name for name in os.listdir('../Spectral-Processing/GCM-OUTPUT/') if os.path.isdir(os.path.join('../Spectral-Processing/GCM-OUTPUT/', name))]
planet_names = ['Taub']

# There are the different sets of opacity and EOS files
# There are somethings that need to be changed in the template inputs file to make this happen
# If you change the underlying data files these might need to be changed
#if any(substring in planet_names[0].upper() for substring in ["GJ1214"]):
#    opacity_files = 'SET_1'
#elif any(substring in planet_names[0].upper() for substring in ["HD209", "HD189"]):
#    opacity_files = 'SET_3'
#else:
#    print("Something is going wrong with how the opacity files are being chosen")
#    exit(0)

opacity_files = 'SET_1'

if opacity_files == 'SET_1':
    cloud_wavelength = 5.0
elif opacity_files == 'SET_2':
    cloud_wavelength = 2.3
elif opacity_files == 'SET_3':
    cloud_wavelength = 2.3
else:
    print("YOU NEED TO SET WHICH OPACITY SET YOU'RE USING")
    exit(0)


# Get the number of GCMS
num_gcms = len(planet_names)

# Get the number of orders of magnitude
num_orders_of_magnitude = grab_input_data.get_input_data('../Spectral-Processing/GCM-OUTPUT/', planet_names[0], "fort.7", "OOM_IN")

planet_radii = []
for planet_name in planet_names:
    planet_radii.append(float(grab_input_data.get_input_data('../Spectral-Processing/GCM-OUTPUT/', planet_name, "fort.7","RADEA")))

# Get the planet name from the input string
planet_name_base = re.split(r"[_|-]", planet_names[0])[0]
planet_name_char_len = len(planet_name_base) + 1

# Auto parse some of these params
# If the file doesn't exist, take a guess and print a warning
column_names = ['lat', 'lon', 'level','alt', 'pres', 'temp', 'u', 'v', 'w']
if os.path.isfile('../Spectral-Processing/PLANET_MODELS/' + planet_names[0] + '.txt'):
    df = pd.read_csv('../Spectral-Processing/PLANET_MODELS/' + planet_names[0] + '.txt', delim_whitespace=True, names=column_names)
    nlat = len(set(df.lat))
    nlon = len(set(df.lon))
    nlev = len(set(df.level))
else:
    print("nlat, nlon, nlev are being set manually in create_all_figures")
    nlat = 48
    nlon = 96
    nlev = 50


# Code to check whether the names are a list
# If not a list, converts it to a list because that's how the plotting functions expect it
if isinstance(planet_names, list):
  pass
else:
    planet_names = [planet_names]

"""
# Plot the broadband phase curves
print ("Plotting the broadband phase curves...")
print ()
print ()
broadband_phase_curves.plot_thermal_phasecurves(planet_names, nlat, nlon, nlev, num_gcms,planet_name_char_len, two_sets_of_planets=False)
broadband_phase_curves.plot_reflected_phasecurves(planet_names, nlat, nlon, nlev, num_gcms,planet_name_char_len, two_sets_of_planets=False)
broadband_phase_curves.plot_reflected_starlight_maps(planet_names, nlat, nlon, nlev, num_gcms, two_sets_of_planets=False)

# Plot the isobaric cloud maps
# If the extra_pressure level is greater than 0, its plots it also
# If it is 0, then it plots the IR photosphere pressure
print ("Plotting the isobaric projections...")
print ()
print ()
cloud_coverage_isobars.plot_cloud_coverage_isobars(planet_names, nlat, nlon, nlev, num_gcms, cloud_wavelength,extra_pressure_level_bar=0)
cloud_coverage_isobars.plot_cloud_coverage_isobars(planet_names, nlat, nlon, nlev, num_gcms, cloud_wavelength,extra_pressure_level_bar=0.001)



# Plot the ptc curves
print ("Plotting the pressure temperature condensation curves...")
print ()
print ()
pressure_temperature_condensation_curves.plot_PTC_curves(planet_names, nlat, nlon, nlev, num_gcms, num_orders_of_magnitude)


# Plot other planet characteristics
aerosol_maps.plot_aerosol_maps(planet_names, nlat, nlon, nlev, num_orders_of_magnitude, cloud_wavelength)
aerosol_profiles.plot_aersol_profiles(planet_names, nlat, nlon, nlev, num_orders_of_magnitude)
wind_maps.plot_wind_maps(planet_names, nlat, nlon, nlev, num_orders_of_magnitude)

# Plotting the emission maps
#emission_maps.plot_emission_maps(planet_names, nlat, nlon, nlev)
"""

# Plot the spectra
print ("Plotting the spectra...")
print ()
print ()
spectra.plot_planet_spectra_blackbody_comparison(planet_names,
                                                 black_body_temperatures=np.linspace(500, 2000, 4),
                                                 num_phases=4)
spectra.plot_star_spectra_test(planet_names)
spectra.plot_filters(planet_names)
spectra.plot_spectra_simple(planet_names, num_phases=24)

# If resolution is set to 0, don't convolve at all
spectra.plot_fp_spectra(planet_names,
                            planet_radii,
                            num_phases=24,
                            transmission_filter_name='None',
                            wav_subset=[2e-6, 3e-6],
                            resolution=500)

spectra.plot_fp_fs_spectra(planet_names,
                            planet_radii,
                            num_phases=24,
                            transmission_filter_name='None',
                            wav_subset=[2e-6, 3e-6],
                            resolution=500)

# If resolution is set to 0, don't convolve at all
spectra.plot_fp_phase_curves(planet_names,
                          planet_name_char_len,
                          num_phases=24,
                          transmission_filter_name='None',
                          wav_subset=[2e-6,3e-6])

spectra.plot_fp_fs_phase_curves(planet_names,
                          planet_name_char_len,
                          planet_radii,
                          num_phases=24,
                          transmission_filter_name='None',
                          wav_subset=[2e-6,3e-6])


