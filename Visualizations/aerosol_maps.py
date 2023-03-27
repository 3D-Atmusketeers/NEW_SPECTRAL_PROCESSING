import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import grab_input_data

def plot_aerosol_maps(planet_names, nlat, nlon, nlev, num_orders_of_magnitude, cloud_wavelength):
    # cloud colormap
    cm_name = 'devon' 
    cm_file = np.loadtxt(f'ScientificColourMaps7/{cm_name}/{cm_name}.txt')
    my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file[0:240])

    def get_optical_depth(file, nlat, nlev):
        df = pd.read_csv(file,
                        delim_whitespace=True, skiprows=0,
                        names=('lat', 'lon', 'level',
                                'alt', 'pres', 'temp', 
                                'u', 'v', 'w',
                                'aero_sw_tau_1', 'sw_asym_1', 'sw_pi0_1',
                                'aero_sw_tau_2', 'sw_asym_2', 'sw_pi0_2',
                                'aero_sw_tau_3', 'sw_asym_3', 'sw_pi0_3',
                                'aero_sw_tau_4', 'sw_asym_4', 'sw_pi0_4',
                                'aero_sw_tau_5', 'sw_asym_5', 'sw_pi0_5',
                                'aero_sw_tau_6', 'sw_asym_6', 'sw_pi0_6',
                                'aero_sw_tau_7', 'sw_asym_7', 'sw_pi0_7',
                                'aero_sw_tau_8', 'sw_asym_8', 'sw_pi0_8',
                                'aero_sw_tau_9', 'sw_asym_9', 'sw_pi0_9',
                                'aero_sw_tau_10', 'sw_asym_10', 'sw_pi0_10',
                                'aero_sw_tau_11', 'sw_asym_11', 'sw_pi0_11',
                                'aero_sw_tau_12', 'sw_asym_12', 'sw_pi0_12',
                                'aero_sw_tau_13', 'sw_asym_13', 'sw_pi0_13',
                                'haze_tau_optical_depth_per_bar', 'haze_asym', 'haze_pi0'))

        df['layer_pressure'] = df['pres'].diff()
        df = df[df.index % nlev != 0]  # Excludes every 250rd row starting from 0

        df["aero_sw_tau_1"] =  df["aero_sw_tau_1"] / df['layer_pressure']
        df["aero_sw_tau_2"] =  df["aero_sw_tau_2"] / df['layer_pressure']
        df["aero_sw_tau_3"] =  df["aero_sw_tau_3"] / df['layer_pressure']
        df["aero_sw_tau_4"] =  df["aero_sw_tau_4"] / df['layer_pressure']
        df["aero_sw_tau_5"] =  df["aero_sw_tau_5"] / df['layer_pressure']
        df["aero_sw_tau_6"] =  df["aero_sw_tau_6"] / df['layer_pressure']
        df["aero_sw_tau_7"] =  df["aero_sw_tau_7"] / df['layer_pressure']
        df["aero_sw_tau_8"] =  df["aero_sw_tau_8"] / df['layer_pressure']
        df["aero_sw_tau_9"] =  df["aero_sw_tau_9"] / df['layer_pressure']
        df["aero_sw_tau_10"] =  df["aero_sw_tau_10"] / df['layer_pressure']
        df["aero_sw_tau_11"] =  df["aero_sw_tau_11"] / df['layer_pressure']
        df["aero_sw_tau_12"] =  df["aero_sw_tau_12"] / df['layer_pressure']
        df["aero_sw_tau_13"] =  df["aero_sw_tau_13"] / df['layer_pressure']


        df = df[(df['lon'] == 0)].reset_index(drop=True)
        df = df.groupby(['lat', 'pres']).mean().reset_index()

        df["total_optical_depth"] = df["aero_sw_tau_1"] + df["aero_sw_tau_2"] + df["aero_sw_tau_3"] + df["aero_sw_tau_4"] + \
                                    df["aero_sw_tau_5"] + df["aero_sw_tau_6"] + df["aero_sw_tau_7"] + df["aero_sw_tau_8"] + \
                                    df["aero_sw_tau_9"] + df["aero_sw_tau_10"] + df["aero_sw_tau_11"] + df["aero_sw_tau_12"] + \
                                    df["aero_sw_tau_13"]

        df[df < 1e-10] = 1e-10
        optical_depth = df["total_optical_depth"].values.reshape(nlat, nlev - 1)
        return optical_depth



    for i in range(len(planet_names)):
        planet_name = planet_names[i]

        molef = grab_input_data.get_input_data('../Spectral-Processing/GCM-OUTPUT/', planet_name, 'fort.7', 'MOLEF')
        hazes_str = grab_input_data.get_input_data('../Spectral-Processing/GCM-OUTPUT/', planet_name, 'fort.7', 'HAZES')[0]

        if hazes_str == 'T':
            hazes = True
        else:
            hazes = False

        if any(i > 1e-20 for i in molef) or hazes:
            plt.clf()
            
            file1 = '../Spectral-Processing/PLANET_MODELS/' + planet_name + '_with_clouds_and_wavelength_dependence.txt'

            fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(9,8), sharex=True, sharey=True)
            plt.subplots_adjust(wspace=0.02, hspace=0.1)

            lats = np.linspace(-90, 90, nlat)
            pressures = np.logspace(2 - num_orders_of_magnitude, 2, nlev - 1)

            optical_depths1 = get_optical_depth(file1, nlat, nlev)
            mp1 = axes.contourf(lats, pressures, np.log10(optical_depths1.T), cmap=my_colors, levels=np.linspace(-2,2,100), extend='both')

            axes.set_yscale("log")
            axes.set_ylim(1e2, 10 ** (2 - num_orders_of_magnitude))

            #cbar = fig.colorbar(mp1, ax=axes.ravel().tolist(), location='top', aspect=50, pad=0.02)
            cbar = fig.colorbar(mp1, aspect=20, pad=0.02)
            cbar.set_label(r'log$_{10}$ Aerosol Optical Depth per bar at ' + str(cloud_wavelength) + ' $\mu$m',
                                                                                                 size=20, labelpad=10)
            cbar.ax.tick_params(labelsize=22)  # set your label size here

            fig.text(0.5, 0.02, r"Latitude (degrees)", size=28, ha='center')
            fig.text(0.0, 0.45, r"Pressure (bar)", size=28, va='center', rotation='vertical')

            plt.savefig('../Figures/Aerosol_Maps_{}.png'.format(planet_name), bbox_inches='tight', dpi=100)
        else:
            pass
    return None
