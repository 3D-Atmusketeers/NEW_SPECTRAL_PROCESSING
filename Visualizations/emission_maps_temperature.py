import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as mcolors
from numpy import ma
from matplotlib import ticker, cm


try:
    from inspect import cleandoc as dedent
except ImportError:
    # Deprecated as of version 3.1. Not quite the same
    # as textwrap.dedent.
    from matplotlib.cbook import dedent


def plot_emission_maps(planet_names, nlat, nlon):
    for planet_name in planet_names:
        file = planet_name

        base = "../Spectral-Processing/FINISHED_SPECTRA/Spec_0_"

        phase = '180.0'

        full_df = pd.read_csv(base + file + "_phase_" + phase + "_inc_0.0000.00_emission_map.dat",
                        names=['tau_index', 'wavelength_m', 'lon', 'lat', 'pressure_pa', 'temp', 'vlos'],
                        sep='\s+')
        
        
        

        wavelengths = list(set(full_df["wavelength_m"]))
        wavelengths = sorted(wavelengths)
        wavelengths = wavelengths[3:4]


        # colormap
        cm_name = 'lajolla'
        cm_file = np.loadtxt(f'ScientificColourMaps8/' + cm_name + '/'+ cm_name +'.txt')
        my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file)

        for i in range(len(wavelengths)):
            plt.clf()
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15,15))
            plt.subplots_adjust(wspace=0.2, hspace=0.02)

            df = full_df[full_df["wavelength_m"] == wavelengths[i]]

            map = Basemap(projection='ortho',lat_0=0,lon_0=180,resolution='l', ax=ax)

            map.drawparallels(np.arange(-90.,90.,30.))
            map.drawmeridians(np.arange(0.,360.,60.))

            # compute native map projection coordinates of lat/lon grid.
            nlon_vis = int(nlon / 2 + 1)
            lons = np.reshape(list(df.lon), (nlat, nlon_vis))
            lats = np.reshape(list(df.lat), (nlat, nlon_vis))

            x, y = map(lons, lats)

            pressure = np.reshape(list(df.pressure_pa), (nlat, nlon_vis))
            temps    = np.reshape(list(df.temp), (nlat, nlon_vis))
            vlos     = np.reshape(list(df.vlos), (nlat, nlon_vis))

            # , locator=ticker.LogLocator()
            # This is in mbar, thats where the 1e2 comes from. 1e5 to Pa, 1e3 to mbar
            #emap = map.contourf(x, y, pressure/1e2,
            #                    levels=20,
            #                    cmap=my_colors)

            emap = map.contourf(x, y, temps,
                                levels=np.linspace(1300, 2000, 100),
                                cmap=my_colors)

            #map.contour(x, y, vlos, levels=[0], colors='white', linewidths=2.0, zorder=3)
            #map.contour(x, y, vlos, levels=[-2000, -1500, -1000, -500], colors='red', alpha=1.0, linewidths=2.0, linestyles='solid', zorder=3)
            #map.contour(x, y, vlos, levels=[500, 1000, 1500, 2000],  alpha=1.0, colors='#2c6fbb', linewidths=2.0,zorder=5)

            wav_str = str(np.round(wavelengths[i] * 1e6, 3))

            cb = fig.colorbar(emap, location='bottom', aspect=15, pad=0.02, fraction=0.06)
            cb.outline.set_linewidth(2)
            cb.set_label(label=r'$\tau$=2/3 Temperature at ' + wav_str + ' $\mu$m', fontsize=28)
            cb.ax.tick_params(labelsize=28) 
            cb.ax.minorticks_on()

            # Set 5 tick labels
            ticks = np.linspace(1300, 2000, 6)  # Generate 5 linearly spaced ticks from 1000 to 2000
            cb.set_ticks(ticks)  # Set these ticks on the colorbar
            cb.set_ticklabels([f"{tick:.0f}" for tick in ticks])  # Optional: format tick labels


            plt.savefig('../Figures/{}_temperature_emission_map_{}_phase_{}.png'.format(file, wav_str, phase), bbox_inches='tight', dpi=200)
