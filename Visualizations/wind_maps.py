import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import rcParams, rc

def plot_wind_maps(planet_names, nlat, nlon, nlev, num_orders_of_magnitude):
    # temp colormap
    cm_name = 'turku'
    cm_name = 'cork'
    cm_file = np.loadtxt(f'ScientificColourMaps7/{cm_name}/{cm_name}.txt')
    cm_file = np.flip(cm_file, axis=0)
    my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file[120:])

    def get_wind_speeds(file):
        df = pd.read_csv(file,
                        delim_whitespace=True, skiprows=0,
                        names=('lat', 'lon', 'level',
                                'alt', 'pres', 'temp', 
                                'u', 'v', 'w',
                                'aero_sw_tau_1', 'sw_asym_1', 'sw_pi0_1'))

        df = df[["lat", "level", "pres", "u"]]
        df = df.groupby(['lat', 'level']).mean().reset_index()
        df = df[["lat", "pres","u"]]

        u_data = df.u.values.reshape(nlat, nlev)
        
        return u_data

    for planet_name in planet_names:
        plt.clf()
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(9,8), sharex=True, sharey=True)
        plt.subplots_adjust(wspace=0.02, hspace=0.1)

        lats = np.linspace(-90, 90, nlat)
        pressures = np.logspace(2 - num_orders_of_magnitude, 2, nlev)

        u_data_1 = get_wind_speeds('../Spectral-Processing/PLANET_MODELS/' + planet_name + '.txt')
        u_data_2 = get_wind_speeds('../Spectral-Processing/PLANET_MODELS/' + planet_name + '.txt')

        mp1 = axes.contourf(lats, pressures, u_data_1.T,  cmap=my_colors, levels=100)
        cb1 = axes.contour(lats, pressures, u_data_1.T, levels=[0], colors='black', linewidths=3, linestyles='dashed')

        #mp2 = axes[1].contourf(lats, pressures, u_data_2.T,  cmap=my_colors, levels=100)
        #cb2 = axes[1].contour(lats, pressures, u_data_2.T, levels=[0], colors='black', linewidths=3, linestyles='dashed')

        axes.set_yscale("log")
        axes.set_ylim(1e2, 10 ** (2 -num_orders_of_magnitude))

        #temp_cbar = fig.colorbar(mp, ax=axes.ravel().tolist(), location='top', aspect=50, pad=0.02)
        temp_cbar = fig.colorbar(mp1, aspect=20, pad=0.02)

        temp_cbar.set_label('E-W Wind Speed (m/s)', size=26)

        fig.text(0.5, 0.05, r"Latitude (degrees)", size=30,  ha='center')
        fig.text(0.0, 0.5, r"Pressure (bar)", size=30,  va='center', rotation='vertical')

        plt.savefig('../Figures/wind_maps_{}.png'.format(planet_name), bbox_inches='tight', dpi=100)
