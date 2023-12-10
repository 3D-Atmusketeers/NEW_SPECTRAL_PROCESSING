#!/usr/bin/env python

# get the data from the fort.7 files
import grab_input_data

### ----- IMPORT LIBRARIES ----- ###
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import ticker
import matplotlib
import matplotlib.colors as colors


def plot_wind_isobars(
        planet_names,
        nlat,
        nlon,
        nlev,
        cloud_wavelength,
        plot_hazes,
        extra_pressure_level_bar):
            
    print('Plotting weird wind isobars')

    # temp colormap
    cm_name = 'cork'
    cm_file = np.loadtxt(f'ScientificColourMaps8/{cm_name}/{cm_name}.txt')
    cm_file = np.flip(cm_file, axis=0)
    temperature_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file)

    # cloud colormap
    cm_name = 'devon'
    cm_file = np.loadtxt(f'ScientificColourMaps8/{cm_name}/{cm_name}.txt')
    cloud_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file[0:240])

    column_names = ['lat', 'lon', 'level',
                    'alt', 'pres', 'temp',
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
                    'haze_tau_optical_depth_per_bar', 'haze_asym', 'haze_pi0']
    nlevel = nlev
    nparams = len(column_names)

    for ind, planet_name in enumerate(planet_names):
        molef = grab_input_data.get_input_data('../Spectral-Processing/GCM-OUTPUT/', planet_name, 'fort.7', 'MOLEF')
        hazes_str = grab_input_data.get_input_data('../Spectral-Processing/GCM-OUTPUT/', planet_name, 'fort.7', 'HAZES')[0]

        if hazes_str == 'T':
            hazes = True
        else:
            hazes = False

        gravity = grab_input_data.get_input_data('../Spectral-Processing/GCM-OUTPUT/', planet_name, 'fort.7', 'GA')
        ir_absorbtion_coefficient = grab_input_data.get_input_data('../Spectral-Processing/GCM-OUTPUT/',
                                                                   planet_name, 'fort.7', 'ABSLW')
        ir_photosphere_pressure_bars = (2. / 3.) * (gravity / ir_absorbtion_coefficient) / 10000
        ir_photosphere_pressure_bars = np.round(ir_photosphere_pressure_bars, 3)

        if (extra_pressure_level_bar == 0):
            pressure_lev_bar = ir_photosphere_pressure_bars
        else:
            pressure_lev_bar = extra_pressure_level_bar

        # Photospheric pressure in bars
        P_phots = [pressure_lev_bar] * len(planet_names)

        plt.clf()
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(12, 9))
        plt.subplots_adjust(wspace=0.01, hspace=0.015)

        P_phot = P_phots[ind]

        df = pd.read_csv(
            '../Spectral-Processing/PLANET_MODELS/' +
            planet_name +
            '_with_clouds_and_wavelength_dependence.txt',
            delim_whitespace=True,
            names=column_names)
        data = df.to_numpy()
        data = data.reshape((nlat, nlon, nlevel, nparams))

        # need to get indices for given pressure level
        pressure_ind = np.zeros((nlat, nlon))

        for i in range(nlat):
            for j in range(nlon):
                k = nlevel - 1

                while data[i][j][k][4] > P_phot:
                    k -= 1
                pressure_ind[i][j] = k

            lons = np.arange(-360, 360, 3.75)
        lats = data[:, :, 0, 0][:, 0]

        ### ----- GET LOCAL TEMP, AEROSOL TAU, etc... ----- ###

        print('Calculating T and aerosol tau...')

        temps = np.zeros((nlat, nlon))
        taus = np.zeros((nlat, nlon))
        EW_vels = np.zeros((nlat, nlon))
        NS_vels = np.zeros((nlat, nlon))
        vert_vels = np.zeros((nlat, nlon))
        pressures = np.zeros((nlat, nlon))
        z = np.zeros((nlat, nlon))

        # empty array to store values
        taus = np.zeros((nlat, nlon))

        for i in range(nlat):
            for j in range(nlon):

                z[i][j] = data[i][j][int(pressure_ind[i][j])][3]
                pressures[i][j] = data[i][j][int(pressure_ind[i][j])][4]
                temps[i][j] = data[i][j][int(pressure_ind[i][j])][5]
                EW_vels[i][j] = data[i][j][int(pressure_ind[i][j])][6]
                NS_vels[i][j] = data[i][j][int(pressure_ind[i][j])][7]
                vert_vels[i][j] = data[i][j][int(pressure_ind[i][j])][8]

                # integrate aerosol optical depth above pressurm level
                k = 0

                while k <= int(pressure_ind[i][j]):
                    taus[i][j] += data[i][j][k][9] + data[i][j][k][12] + data[i][j][k][15] + data[i][j][k][18] + \
                                  data[i][j][k][21] + data[i][j][k][24] + data[i][j][k][27] + data[i][j][k][30] + \
                                  data[i][j][k][33] + data[i][j][k][36] + data[i][j][k][39] + data[i][j][k][42] + \
                                  data[i][j][k][45]


                    # The hazes are optical depth per bar
                    if plot_hazes:
                        taus[i][j] += data[i][j][k][48] * (data[i][j][k+1][4] - data[i][j][k][4])
                    k += 1

        # stream plot
        test_x = np.linspace(min(lons), max(lons), len(lons))
        test_y = np.linspace(min(lats), max(lats), len(lats))


        
        if (pressure_lev_bar < 1e-2):
            my_levels = np.arange(700, 1800, 100)
        elif 'hd189' in planet_name.lower():
            my_levels = np.arange(600, 1600, 100)
        else:
            my_levels = np.arange(1000, 1800, 100)

        matplotlib_version = matplotlib.__version__.split('.')
        if float(matplotlib_version[1]) < 2:
            my_norm = colors.DivergingNorm(0)
        else:
            my_norm = colors.TwoSlopeNorm(0)

        temp_map = axes.contourf(lons,
                                 lats, np.concatenate([EW_vels, EW_vels], axis=1),
                                 cmap=temperature_colors,
                                 norm=my_norm)
        temp_cbar = fig.colorbar(
        temp_map,
        aspect=30,
        pad=0.015,
        orientation='horizontal')
            
        temp_cbar.set_label('E-W Wind Speed (m/s)', fontsize=26)
        
        axes.streamplot(test_x, test_y,
                               np.concatenate([EW_vels, EW_vels], axis=1),
                               np.concatenate([NS_vels, NS_vels], axis=1),
                               linewidth=1.2, density=([1, 1.5]), color='#d8dcd6', zorder=1)

        # format axes
        axes.set_xlim([-180, 180])
        axes.set_ylim([-87, 87])
        axes.spines['left'].set_position('zero')
        axes.spines['bottom'].set_position('zero')
        axes.spines['left'].set_color('none')
        axes.spines['right'].set_color('none')
        axes.spines['top'].set_color('none')
        axes.spines['bottom'].set_color('none')
        axes.tick_params(
            axis='both',
            which='both',
            top=False,
            bottom=False,
            right=False,
            left=False,
            labelcolor='w',
            labelsize=25,
            pad=-12)

        # axes.yaxis.get_major_ticks()[2].label1.set_visible(False)
        axes.set_yticks([-75, -50, -25, 25, 50, 75])
        axes.grid(color='w', alpha=0.5, ls=':')

        axes.patch.set_edgecolor('black')
        axes.patch.set_linewidth('2')
        axes.grid(color='w', alpha=0.5, ls=':')


        plt.savefig('../Figures/Wind_Isobars_{}_bar_{}.png'.format(
            P_phots[0], planet_name), bbox_inches='tight', dpi=250)
