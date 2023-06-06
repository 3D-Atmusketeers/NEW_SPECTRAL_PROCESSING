import pandas as pd
import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import grab_input_data


def add_vlos(file, r_p, omega, inc):
    df = pd.read_csv(file + '.txt', delim_whitespace=True, names=('lat', 'lon', 'level',
                                         'alt', 'pressure', 'temp',
                                         'u', 'v', 'w',
                                         'aero_tau_pre_qext_1', 'sw_asym_1', 'sw_pi0_1',
                                         'aero_tau_pre_qext_2', 'sw_asym_2', 'sw_pi0_2',
                                         'aero_tau_pre_qext_3', 'sw_asym_3', 'sw_pi0_3',
                                         'aero_tau_pre_qext_4', 'sw_asym_4', 'sw_pi0_4',
                                         'aero_tau_pre_qext_5', 'sw_asym_5', 'sw_pi0_5',
                                         'aero_tau_pre_qext_6', 'sw_asym_6', 'sw_pi0_6',
                                         'aero_tau_pre_qext_7', 'sw_asym_7', 'sw_pi0_7',
                                         'aero_tau_pre_qext_8', 'sw_asym_8', 'sw_pi0_8',
                                         'aero_tau_pre_qext_9', 'sw_asym_9', 'sw_pi0_9',
                                         'aero_tau_pre_qext_10', 'sw_asym_10', 'sw_pi0_10',
                                         'aero_tau_pre_qext_11', 'sw_asym_11', 'sw_pi0_11',
                                         'aero_tau_pre_qext_12', 'sw_asym_12', 'sw_pi0_12',
                                         'aero_tau_pre_qext_13', 'sw_asym_13', 'sw_pi0_13',
                                         'haze_tau_optical_depth_per_bar', 'haze_asym', 'haze_pi0',
                                         'incident_frac'))
    dtr = 0.0174533 # Degrees to radians
    df['vlos'] =    df.u * np.sin(df.lon * dtr) \
                  + df.v * np.cos(df.lon * dtr) * np.sin(df.lat * dtr) \
                  + np.cos(inc) * (omega * (r_p + df.alt) * np.sin(df.lon * dtr) * np.cos(df.lat * dtr))
    return df



def plot_observer_projection(planet_names, nlat, nlon,
                             planet_radii, pressure_in_mbar=10):
    
    for i, planet_name in enumerate(planet_names):
        planet_radius = planet_radii[i]

        # Colormap
        cm_name = 'lajolla'
        cm_file = np.loadtxt('ScientificColourMaps7/lajolla/lajolla.txt')
        my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file[::-1])

        # Create figure and adjust subplots
        fig, ax = plt.subplots(figsize=(23, 24))
        plt.subplots_adjust(wspace=-0.2, hspace=0.02)

        # File path
        file = '../Spectral-Processing/Spectra/DATA/init_' + planet_name + '_phase_0.0_inc_0.0'

        rotation_rate = grab_input_data.get_input_data('../Spectral-Processing/GCM-OUTPUT/', planet_name, 'fort.7', 'WW')

        print(rotation_rate)

        # Parameters
        pressure_min = pressure_in_mbar
        r_p = planet_radius
        omega = rotation_rate
        inc = 0.0

        # Calculate vlos and other variables
        df = add_vlos(file, r_p, omega, inc)
        all_lats = sorted(list(set(df.lat)))
        all_lons = sorted(list(set(df.lon)))

        lats, lons, temps, us, vs, vloss, taus = [], [], [], [], [], [], []

        for lat in all_lats:
            for lon in all_lons:
                temp_df = df[(df.lat == lat) & (df.lon == lon)].reset_index()
                
                k = np.argmax(temp_df.pressure * 1e-2 >= pressure_min)
                
                lats.append(lat)
                lons.append(lon)
                temps.append(temp_df.temp[k])
                vloss.append(temp_df.vlos[k])
                us.append(temp_df.u[k])
                vs.append(temp_df.v[k])
                taus.append(0)

        # Reshape variables
        new_lats = np.reshape(lats, (nlat, nlon))
        new_lons = np.reshape(lons, (nlat, nlon))
        final_lats = new_lats * np.pi / 180
        final_lons = new_lons * np.pi / 180
        final_taus = np.reshape(taus, (nlat, nlon))
        final_temps = np.reshape(temps, (nlat, nlon))
        final_us = np.reshape(us, (nlat, nlon))
        final_vs = np.reshape(vs, (nlat, nlon))
        final_vloss = np.reshape(vloss, (nlat, nlon))

        # Create basemap
        map = Basemap(projection='ortho', lat_0=0, lon_0=180, resolution='l')

        # Compute map projection coordinates
        x, y = map(final_lons * 180. / np.pi, final_lats * 180. / np.pi)

        # Plot contours and quiver
        cs1 = map.contourf(x, y, final_temps, 100, cmap=my_colors)
        map.contour(x, y, final_vloss, levels=[0], colors='white', linewidths=3, zorder=3)
        map.contour(x, y, final_vloss, levels=[-2000, -1500, -1000, -500], colors='red', alpha=0.8, linewidths=3, linestyles='solid', zorder=3)
        map.contour(x, y, final_vloss, levels=[500, 1000, 1500, 2000],  alpha=0.8, colors='navy', linewidths=3,zorder=5)

        # Graph the wind vectors
        #map.quiver(x[::2,::2], y[::2,::2], final_us[::2,::2], final_vs[::2,::2], color='black', zorder=6) 

        cbar1 = fig.colorbar(cs1, location='top', aspect=40, pad=0.02)
        cbar1.outline.set_linewidth(2)
        cbar1.set_label(label='Temperature (K)',weight='bold', fontsize=25)
        cbar1.ax.tick_params(labelsize=20) 

        fig.savefig('../Figures/observer_projection_{}.png'.format(planet_name), bbox_inches='tight', dpi=100)
