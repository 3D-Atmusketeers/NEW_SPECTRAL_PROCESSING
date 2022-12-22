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

import aerosol_maps
import aerosol_profiles
import wind_maps
import emission_maps

plt.style.use('style.txt')

planet_names = [name for name in os.listdir('../Spectral-Processing/GCM-OUTPUT/') if os.path.isdir(os.path.join('../Spectral-Processing/GCM-OUTPUT/', name))]
planet_names = ['HD189-DOGRAY-NUC-CLOUDS-COMPACT']

# Get the planet name from the input string
planet_name_base = re.split(r"[_|-]", planet_names[0])[0]

# Get the name of the planet name string, add 1 to account for the - or _
planet_name_char_len = len(planet_name_base) + 1 

column_names = ['lat', 'lon', 'level','alt', 'pres', 'temp', 'u', 'v', 'w']
df = pd.read_csv('../Spectral-Processing/PLANET_MODELS/' + planet_names[0] + '.txt', delim_whitespace=True, names=column_names)

nlat = len(set(df.lat))
nlon = len(set(df.lon))
nlev = len(set(df.level))

num_gcms = len(planet_names)

# The number of orders of magnitude in the GCM, should really be read in from the fort file
num_orders_of_magnitude = 7

# Code to check whether the names are a list
# If not a list, converts it to a list because that's how the plotting functions expect it
if isinstance(planet_names, list):
  pass
else:
    planet_names = [planet_names]

#print ("Plotting the broadband phase curves...")
#print ()
#print ()
# Plot the broadband phase curves
#broadband_phase_curves.plot_thermal_phasecurves(planet_names, nlat, nlon, nlev, num_gcms,planet_name_char_len)
#broadband_phase_curves.plot_reflected_phasecurves(planet_names, nlat, nlon, nlev, num_gcms,planet_name_char_len)
#broadband_phase_curves.plot_reflected_starlight_maps(planet_names, nlat, nlon, nlev, num_gcms)

#print ("Plotting the isobaric projections...")
#print ()
#print ()
# Plot the isobaric cloud maps
# If the extra_pressure level is greater than 0, its plots it also
# If it is 0, then it plots the IR photosphere pressure
#cloud_coverage_isobars.plot_cloud_coverage_isobars(planet_names, nlat, nlon, nlev, num_gcms, extra_pressure_level_bar=0)
#cloud_coverage_isobars.plot_cloud_coverage_isobars(planet_names, nlat, nlon, nlev, num_gcms, extra_pressure_level_bar=0.1)

#print ("Plotting the pressure temperature condensation curves...")
#print ()
#print ()

# Plot the ptc curves
#pressure_temperature_condensation_curves.plot_PTC_curves(planet_names, nlat, nlon, nlev, num_gcms, nucleation_lim=True)

#print ("Plotting the spectra...")
#print ()
#print ()

# Plot the spectra
spectra.plot_planet_spectra_blackbody_comparison(planet_names, black_body_temperatures=np.linspace(500, 1000, 6), num_phases=1)
#spectra.plot_star_spectra_test(planet_names)
#spectra.plot_filters(planet_names)
#spectra.plot_spectra_simple(planet_names, num_phases=24)
#spectra.plot_spectra_phases(planet_names, num_phases=24, transmission_filter_name='MIRI', planet_only_bool=False)
#spectra.plot_phase_curves(planet_names, planet_name_char_len, num_phases=24, transmission_filter_name='MIRI')

# Plot other stuff
#aerosol_maps.plot_aerosol_maps(planet_names, nlat, nlon, nlev)
#aerosol_profiles.plot_aersol_profiles(planet_names, nlat, nlon, nlev)
#wind_maps.plot_wind_maps(planet_names, nlat, nlon, nlev, num_orders_of_magnitude)
#emission_maps.plot_emission_maps(planet_names, nlat, nlon, nlev)