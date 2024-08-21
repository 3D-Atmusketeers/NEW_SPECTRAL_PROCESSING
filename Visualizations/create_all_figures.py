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
import aerosol_coverage_isobars
import pressure_temperature_condensation_curves
import spectra

import grab_input_data

import aerosol_maps
import aerosol_profiles
import wind_maps
import wind_isobars
import emission_maps_pressure
import emission_maps_temperature
import basemap_hemispheric_projections

import cross_correlation


import scienceplots
plt.style.use('science.mplstyle')
#plt.style.use('astrophysics.mplstyle')
#plt.style.use(['science','nature'])

# Figure out what planets!
planet_names = [name for name in os.listdir('../Spectral-Processing/GCM-OUTPUT/') if os.path.isdir(os.path.join('../Spectral-Processing/GCM-OUTPUT/', name))]

planet_names = [
                "GJ1214b-none-0clouds-1met",
                "GJ1214b-soot-0clouds-1met",
                "GJ1214b-tholin-0clouds-1met",
                "GJ1214b-soot_2xpi0-0clouds-1met",

                "GJ1214b-none-0clouds-30met",
                "GJ1214b-soot-0clouds-30met",
                "GJ1214b-tholin-0clouds-30met",
                "GJ1214b-soot_2xpi0-0clouds-30met",

                "GJ1214b-none-0clouds-100met",
                "GJ1214b-soot-0clouds-100met",
                "GJ1214b-tholin-0clouds-100met",
                "GJ1214b-soot_2xpi0-0clouds-100met",

                "GJ1214b-none-25clouds-1met",
                "GJ1214b-soot-25clouds-1met",
                "GJ1214b-tholin-25clouds-1met",
                "GJ1214b-soot_2xpi0-25clouds-1met",

                "GJ1214b-none-25clouds-30met",
                "GJ1214b-soot-25clouds-30met",
                "GJ1214b-tholin-25clouds-30met",
                "GJ1214b-soot_2xpi0-25clouds-30met",

                "GJ1214b-none-25clouds-100met",
                "GJ1214b-soot-25clouds-100met",
                "GJ1214b-tholin-25clouds-100met",
                "GJ1214b-soot_2xpi0-25clouds-100met",

                "GJ1214b-none-50clouds-1met",
                "GJ1214b-soot-50clouds-1met",
                "GJ1214b-tholin-50clouds-1met",
                "GJ1214b-soot_2xpi0-50clouds-1met",

                "GJ1214b-none-50clouds-30met",
                "GJ1214b-soot-50clouds-30met",
                "GJ1214b-tholin-50clouds-30met",
                "GJ1214b-soot_2xpi0-50clouds-30met",

                "GJ1214b-none-50clouds-100met",
                "GJ1214b-soot-50clouds-100met",
                "GJ1214b-tholin-50clouds-100met",
                "GJ1214b-soot_2xpi0-50clouds-100met",
                ]


planet_names = [
                "GJ1214b-none-0clouds-1met",
                "GJ1214b-none-0clouds-30met",
                "GJ1214b-none-0clouds-100met",
                ]



#inclination needs to be a string in radians (sorry)
inclination = '0.00'

# Set the cloud wavelength to use in microns
cloud_wavelength = 5.0

# Get the number of GCMS
num_gcms = len(planet_names)

# Get the planet name from the input string
planet_name_base = re.split(r"[_|-]", planet_names[0])[0]
planet_name_char_len = len(planet_name_base) + 1


planet_radii = []
for planet_name in planet_names:
    if os.path.isfile('../Spectral-Processing/GCM-OUTPUT/' + planet_name + '/Planet_Run/fort.7'):
        planet_radii.append(float(grab_input_data.get_input_data('../Spectral-Processing/GCM-OUTPUT/',
                                                            planet_name, "fort.7","RADEA")))
        print("Planet Radii: ", planet_radii[0])
    else:
        if "hd189" in planet_name.lower():
            planet_radii.append(0.8e+08) # HD 189
        elif "hd209" in planet_name.lower():
            planet_radii.append(0.99e+08) # HD 209
        elif "gj1214" in planet_name.lower():
            planet_radii.append(17469282.0)
        elif "wasp-77" in planet_name.lower():
            planet_radii.append(8.603e7)
        elif "wasp-121" in planet_name.lower():
            planet_radii.append(1.304e8)
        else:
            print("WARNING! fort.7 file not found for " + planet_name)
            exit(0)

    
# Auto parse some of these params
# If the file doesn't exist, take a guess and print a warning
print("nlat, nlon, nlev are being set manually in create_all_figures")
column_names = ['lat', 'lon', 'level','alt', 'pres', 'temp', 'u', 'v', 'w']
nlat = 48
nlon = 96
nlev = 50
num_orders_of_magnitude = 7


# Code to check whether the names are a list
# If not a list, converts it to a list because that's how the plotting functions expect it
if isinstance(planet_names, list):
  pass
else:
    planet_names = [planet_names]


# Plot the broadband phase curves
#broadband_phase_curves.plot_reflected_phasecurves(planet_names, nlon, two_sets_of_planets=False)
#broadband_phase_curves.plot_thermal_phasecurves(planet_names, nlon, two_sets_of_planets=False)
#broadband_phase_curves.plot_reflected_starlight_maps(planet_names)


# Plot the isobaric cloud maps
# If the extra_pressure level is greater than 0, its plots it also
# If it is 0, then it plots the IR photosphere pressure
#aerosol_coverage_isobars.plot_aerosol_coverage_isobars(planet_names, nlat, nlon, nlev, cloud_wavelength, plot_hazes=False, extra_pressure_level_bar=0)
#aerosol_coverage_isobars.plot_aerosol_coverage_isobars(planet_names, nlat, nlon, nlev, cloud_wavelength, plot_hazes = True, extra_pressure_level_bar = 1e-3)


# Plot the ptc curves
pressure_temperature_condensation_curves.plot_PTC_curves(planet_names, nlat, nlon, nlev, num_orders_of_magnitude)


# Plot other planet characteristics
#aerosol_maps.plot_aerosol_maps(planet_names, nlat, nlon, nlev, num_orders_of_magnitude, cloud_wavelength)
#aerosol_profiles.plot_aersol_profiles(planet_names, nlat, nlon, nlev, num_orders_of_magnitude)
#wind_maps.plot_wind_maps(planet_names, nlat, nlon, nlev, num_orders_of_magnitude)
#wind_isobars.plot_wind_isobars(planet_names, nlat, nlon, nlev, cloud_wavelength, plot_hazes=False, extra_pressure_level_bar=0.01)


# Plotting the emission maps
#emission_maps_pressure.plot_emission_maps(planet_names, nlat, nlon)
#emission_maps_temperature.plot_emission_maps(planet_names, nlat, nlon)
#basemap_hemispheric_projections.plot_observer_projection(planet_names, nlat, nlon,planet_radii, pressure_in_mbar=10)
#cross_correlation.plot_cross_correlations(planet_names, num_phases=24)

# Plot the spectra
#planet_name = planet_names[0]

######## ######## ######## ######## ######## ######## ######## ######## ######## #
########      These are general plotting functions for testing            ########
# ######## ######## ######## ######## ######## ######## ######## ######## ######## 

#spectra.plot_blackbody_phase_curve(planet_name,planet_radii,num_phases=4,transmission_filter_name='MIRI',wav_subset=[5e-6, 12e-6],resolution=100,temp=600)

#spectra.plot_planet_spectra_blackbody_comparison_hz(planet_names,black_body_temperatures=[553], num_phases=4, inclination)
#spectra.plot_planet_spectra_blackbody_comparison_microns(planet_names,black_body_temperatures=[553],num_phases=2, inclination)
#spectra.plot_star_spectra_test(planet_names)
#spectra.plot_spectra_simple(planet_names, num_phases=4, inclination)


# The options are "SPITZER_3_6, SPITZER_4_5, MIRI"
for filter_name in ['None']:
    if filter_name != 'None':
        print()
        print(filter_name)

        #spectra.plot_filters(planet_names, transmission_filter_name=filter_name, inclination)

    #spectra.plot_fp_spectra(planet_names,
    #                        num_phases=2,
    #                        transmission_filter_name=filter_name,
    #                        wav_subset=[0, 1],
    #                        resolution=100,
    #                        INC_STRING=inclination)

    #spectra.plot_fp_fs_spectra(planet_names,
    #                            planet_radii,
    #                            num_phases=24,
    #                            transmission_filter_name=filter_name,
    #                            wav_subset=[5e-6, 12e-6],
    #                            resolution=100,
    #                            INC_STRING=inclination)

    #spectra.plot_dayside(planet_names,
    #                            planet_radii,
    #                            num_phases=2,
    #                            transmission_filter_name=filter_name,
    #                            wav_subset=[0,12e-6],
    #                            resolution=100,
    #                            INC_STRING=inclination)

    spectra.plot_fp_phase_curves(planet_names,
                            planet_name_char_len,
                            num_phases=24,
                            transmission_filter_name=filter_name,
                            wav_subset=[0, 100],
                            INC_STRING=inclination)

    #spectra.plot_fp_fs_phase_curves(planet_names,
    #                        planet_name_char_len,
    #                        planet_radii,
    #                        num_phases=24,
    #                        transmission_filter_name=filter_name,
    #                        wav_subset=[5e-6, 12.0e-6],
    #                        INC_STRING=inclination)

# Notes
# if the transmission filter name is set to 'None', then it will plot the full spectrum
# unless you specify a wavelength subset
# the options for filters are MIRI, SPITZER_3_6, SPITZER_4_5
# If resolution is set to 0, don't convolve at all

