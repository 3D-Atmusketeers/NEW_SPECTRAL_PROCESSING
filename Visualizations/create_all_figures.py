#!/usr/bin/python
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt

# These are the plotting codes that Isaac wrote
# If they're broken, venmo isaac-malsky
# I wont' fix them, but I'll appreciate the money
import broadband_phase_curves
import cloud_coverage_isobars
import pressure_temperature_condensation_curves
import spectra

plt.style.use('style.txt')

base = '../Spectral-Processing/PLANET_MODELS/'
#planet_names = [name for name in os.listdir(base) if os.path.isdir(os.path.join(base, name))]
planet_names = ["GJ1214b-SOOT-HAZES-1X"]

column_names = ['lat', 'lon', 'level',
               'alt', 'pres', 'temp', 
               'u', 'v', 'w']

df = pd.read_csv(base + planet_names[0] + '.txt', delim_whitespace=True, names=column_names)

nlat = len(set(df.lat))
nlon = len(set(df.lon))
nlev = len(set(df.level))
num_gcms = len(planet_names)


print ("Plotting the broadband phase curves...")
print ()
print ()
# Plot the broadband phase curves
broadband_phase_curves.plot_phasecurves(planet_names, nlat, nlon, nlev, num_gcms)
broadband_phase_curves.plot_reflected_starlight_maps(planet_names, nlat, nlon, nlev, num_gcms)

print ("Plotting the isobaric projections...")
print ()
print ()
# Plot the isobaric cloud maps
pressure_lev_bar = 1e-2
cloud_coverage_isobars.plot_cloud_coverage_isobars(planet_names, nlat, nlon, nlev, num_gcms, pressure_lev_bar)


print ("Plotting the pressure temperature condensation curves...")
print ()
print ()
# Plot the ptc curves
pressure_temperature_condensation_curves.plot_PTC_curves(planet_names, nlat, nlon, nlev, num_gcms, nucleation_lim=True)


print ("Plotting the spectra...")
print ()
print ()

# Plot the spectra
spectra.plot_planet_spectra_test(planet_names)
spectra.plot_star_spectra_test(planet_names)
spectra.plot_filters(planet_names)
spectra.plot_spectra_phases(planet_names)
spectra.plot_phase_curves(planet_names)