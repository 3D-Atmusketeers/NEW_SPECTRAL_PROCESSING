import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as mcolors

try:
    from inspect import cleandoc as dedent
except ImportError:
    # Deprecated as of version 3.1. Not quite the same
    # as textwrap.dedent.
    from matplotlib.cbook import dedent

def plot_emission_maps(planet_names, nlat, nlon, nlev):
    for planet_name in planet_names:
        file = planet_name

        base = "../Spectral-Processing/FINISHED_SPECTRA/Spec_0_"

        full_df = pd.read_csv(base + file + "_phase_0.0_inc_0.00.00.0000.00_emission_map.dat",
                        names=['tau_index', 'wavelength_m', 'lon', 'lat', 'pressure_pa'],
                        delim_whitespace=True)

        wavelengths = list(set(full_df["wavelength_m"]))

        # colormap
        cm_name = 'lapaz'
        cm_file = np.loadtxt(f'ScientificColourMaps7/' + cm_name + '/'+ cm_name +'.txt')
        my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file)

        min_val = (int(np.amin(full_df.pressure_pa/1e3))) / 100
        max_val = (int(np.amax(full_df.pressure_pa/1e3))) / 100

        print (min_val, max_val)


        for i in range(len(wavelengths)):
            plt.clf()
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,10))
            plt.subplots_adjust(wspace=0.2, hspace=0.02)
            
            df = full_df[full_df["wavelength_m"] == wavelengths[i]]
            
            map = Basemap(projection='ortho',lat_0=0,lon_0=180,resolution='l', ax=ax)
            
            map.drawparallels(np.arange(-90.,90.,30.))
            map.drawmeridians(np.arange(0.,360.,60.))

            # compute native map projection coordinates of lat/lon grid.
            lons = np.reshape(list(df.lon), (48,49))
            lats = np.reshape(list(df.lat), (48,49))

            x, y = map(lons, lats)

            pressure = np.reshape(list(df.pressure_pa), (48, 49))

            emap = map.contourf(x, y, pressure/1e5, levels=200, cmap=my_colors)
            #map.contour(x, y, pressure/1e5, levels=[0.2, 0.25, 0.30, 0.35, 0.40], cmap=my_colors)

            wav_str = str(np.round(wavelengths[i] * 1e6, 3))
            cb = map.colorbar(emap, ticks=np.linspace(min_val, max_val, 7),
                            location='bottom', label=r'$\tau$=2/3 Pressure at ' + wav_str + ' $\mu$m (bar)')
            plt.savefig('../Figures/{}_emission_map_{}.png'.format(file, wav_str), bbox_inches='tight', dpi=200)