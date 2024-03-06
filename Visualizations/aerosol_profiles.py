import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import grab_input_data
import matplotlib.pylab as pl
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoMinorLocator

def plot_aersol_profiles(planet_names, nlat, nlon, nlev, num_orders_of_magnitude):

    colors=['#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#ffffff']
    colors = pl.cm.nipy_spectral(np.linspace(0, 1, 13))


    kcl_100x_condensation_curve = [653.809,657.406,662.437,668.463,675.364,682.441,
                          689.268,695.987,703.105,710.176,717.652,725.173,
                          733.103,741.058,749.31,757.571,766.093,774.665,
                          783.505,792.484,801.78,811.299,821.075,831.095,
                          841.321,851.81,862.502,873.611,884.895,896.547,
                          908.319,920.64,933.102,946.178,959.34,973.242,
                          987.244,1002.033,1016.74,1032.244,1047.697,1064.283,
                          1080.85,1098.632,1116.36,1137.157,1158.206,1178.101,
                          1194.859,1208.199]


    for planet_name in planet_names:
        gravity = grab_input_data.get_input_data('../Spectral-Processing/GCM-OUTPUT/', planet_name,'fort.7' ,'GA')

        ir_absorbtion_coefficient = grab_input_data.get_input_data('../Spectral-Processing/GCM-OUTPUT/', planet_name,'fort.7' ,'ABSLW')
        ir_photosphere_pressure_bars = (2./3.) * (gravity/ir_absorbtion_coefficient) / 10000
        ir_photosphere_pressure_bars = np.round(ir_photosphere_pressure_bars, 3)

        vis_absorbtion_coefficient = grab_input_data.get_input_data('../Spectral-Processing/GCM-OUTPUT/', planet_name,'fort.7' ,'ABSSW')
        vis_photosphere_pressure_bars = (2./3.) * (gravity/vis_absorbtion_coefficient) / 10000
        vis_photosphere_pressure_bars = np.round(vis_photosphere_pressure_bars, 3)

        i = 0
        file = '../Spectral-Processing/PLANET_MODELS/' + planet_name +  '_with_clouds_and_wavelength_dependence.txt'
        df1 = pd.read_csv(file,
                        delim_whitespace=True, skiprows=0,
                        names=('lat', 'lon', 'level',
                                            'alt','pres','temp',
                                            'u', 'v', 'w',
                                            'aero_tau_1', 'sw_asym_1', 'sw_pi0_1',
                                            'aero_tau_2', 'sw_asym_2', 'sw_pi0_2',
                                            'aero_tau_3', 'sw_asym_3', 'sw_pi0_3',
                                            'aero_tau_4', 'sw_asym_4', 'sw_pi0_4',
                                            'aero_tau_5', 'sw_asym_5', 'sw_pi0_5',
                                            'aero_tau_6', 'sw_asym_6', 'sw_pi0_6',
                                            'aero_tau_7', 'sw_asym_7', 'sw_pi0_7',
                                            'aero_tau_8', 'sw_asym_8', 'sw_pi0_8',
                                            'aero_tau_9', 'sw_asym_9', 'sw_pi0_9',
                                            'aero_tau_10', 'sw_asym_10', 'sw_pi0_10',
                                            'aero_tau_11', 'sw_asym_11', 'sw_pi0_11',
                                            'aero_tau_12', 'sw_asym_12', 'sw_pi0_12',
                                            'aero_tau_13', 'sw_asym_13', 'sw_pi0_13',
                                            'haze_tau_optical_depth', 'haze_asym', 'haze_pi0'))

        df2 = pd.read_csv(file,
                        delim_whitespace=True, skiprows=0,
                        names=('lat', 'lon', 'level',
                                'alt','pres','temp',
                                'u', 'v', 'w',
                                'aero_tau_1', 'sw_asym_1', 'sw_pi0_1',
                                'aero_tau_2', 'sw_asym_2', 'sw_pi0_2',
                                'aero_tau_3', 'sw_asym_3', 'sw_pi0_3',
                                'aero_tau_4', 'sw_asym_4', 'sw_pi0_4',
                                'aero_tau_5', 'sw_asym_5', 'sw_pi0_5',
                                'aero_tau_6', 'sw_asym_6', 'sw_pi0_6',
                                'aero_tau_7', 'sw_asym_7', 'sw_pi0_7',
                                'aero_tau_8', 'sw_asym_8', 'sw_pi0_8',
                                'aero_tau_9', 'sw_asym_9', 'sw_pi0_9',
                                'aero_tau_10', 'sw_asym_10', 'sw_pi0_10',
                                'aero_tau_11', 'sw_asym_11', 'sw_pi0_11',
                                'aero_tau_12', 'sw_asym_12', 'sw_pi0_12',
                                'aero_tau_13', 'sw_asym_13', 'sw_pi0_13',
                                'haze_tau_optical_depth', 'haze_asym', 'haze_pi0'))

        df1['layer_pressure'] = df1['pres'].diff()
        df2['layer_pressure'] = df2['pres'].diff()

        layer_pressures = list(df1['layer_pressure'].head(50))
        layer_pressures[0] = 10.0 ** (np.log10(layer_pressures[1]) - (np.log10(layer_pressures[2]) - np.log10(layer_pressures[1])))
        layer_pressures = np.asarray(layer_pressures)



        df1 = df1[(df1['lat'] == 1.8556)  & (df1['lon'] == 0.0)].reset_index(drop=True)
        df2 = df2[(df2['lat'] == 1.8556)  & (df2['lon'] == 180.0)].reset_index(drop=True)

        plt.clf()

        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(30,8), sharex=True)
        plt.subplots_adjust(wspace=0.2, hspace=0.15)

        optical_depths_1 = df1.aero_tau_1 + df1.aero_tau_2 + df1.aero_tau_3 + df1.aero_tau_4 + \
                        df1.aero_tau_5 + df1.aero_tau_6 + df1.aero_tau_7 + df1.aero_tau_8 + \
                        df1.aero_tau_9 + df1.aero_tau_10 + df1.aero_tau_11 + df1.aero_tau_12 + \
                        df1.aero_tau_13

        optical_depths_2 = df2.aero_tau_1 + df2.aero_tau_2 + df2.aero_tau_3 + df2.aero_tau_4 + \
                        df2.aero_tau_5 + df2.aero_tau_6 + df2.aero_tau_7 + df2.aero_tau_8 + \
                        df2.aero_tau_9 + df2.aero_tau_10 + df2.aero_tau_11 + df2.aero_tau_12 + \
                        df2.aero_tau_13

        haze_optical_depths_1 = df1.haze_tau_optical_depth
        haze_optical_depths_2 = df2.haze_tau_optical_depth

        cum_optical_depths_1 = np.cumsum(list(optical_depths_1))
        cum_optical_depths_2 = np.cumsum(list(optical_depths_2))

        cum_haze_optical_depths_1 = np.cumsum(list(haze_optical_depths_1))
        cum_haze_optical_depths_2 = np.cumsum(list(haze_optical_depths_2))


        ax[0].plot(df1.pres, haze_optical_depths_1 / layer_pressures, color='red', linewidth=2, linestyle='dashed', label=r'Haze,substellar profile')
        ax[0].plot(df2.pres, haze_optical_depths_2 / layer_pressures, color='black', linewidth=2, linestyle='dashed',label=r'Haze, antistellar profile')

        ax[0].plot(df1.pres, df1.aero_tau_1 / layer_pressures, color='red', linewidth=2, label='KCl, substellar profile')
        ax[0].plot(df2.pres, df2.aero_tau_1 / layer_pressures, color='black', linewidth=2, label='KCl, antistellar profile')

        ax[1].plot(df1.pres, cum_optical_depths_1, color='red', linewidth=2)
        ax[1].plot(df2.pres, cum_optical_depths_2, color='black', linewidth=2)

        ax[1].plot(df2.pres, cum_haze_optical_depths_1, color='red', linewidth=2,linestyle='dashed')
        ax[1].plot(df2.pres, cum_haze_optical_depths_2, color='black', linewidth=2,linestyle='dashed')

        """
        np.savetxt('y_values_haze_sub.txt', haze_optical_depths_1 / layer_pressures, delimiter=',')
        np.savetxt('y_values_haze_anti.txt', haze_optical_depths_2 / layer_pressures, delimiter=',')
        np.savetxt('y_values_kcl_sub.txt', df1.aero_tau_1 / layer_pressures, delimiter=',')
        np.savetxt('y_values_kcl_anti.txt', df2.aero_tau_1 / layer_pressures, delimiter=',')
        np.savetxt('y_values_cum_optical_depths_1.txt', cum_optical_depths_1, delimiter=',')
        np.savetxt('y_values_cum_optical_depths_2.txt', cum_optical_depths_2, delimiter=',')
        np.savetxt('y_values_cum_haze_optical_depths_1.txt', cum_haze_optical_depths_1, delimiter=',')
        np.savetxt('y_values_cum_haze_optical_depths_2.txt', cum_haze_optical_depths_2, delimiter=',')
        """

        ax[2].plot(df1.pres, df1.temp, color='red', linewidth=2, label='Substellar profile')
        ax[2].plot(df2.pres, df2.temp, color='black', linewidth=2, label='Antistellar profile')

        ax[2].plot(df1.pres, kcl_100x_condensation_curve, color='#03719c',
                   linestyle='dashdot',label='KCl Condensation Curve', linewidth=3)

        ax[0].set_xscale('log')
        ax[1].set_xscale('log')
        ax[2].set_xscale('log')

        ax[0].set_yscale('log')
        ax[1].set_yscale('log')
        ax[2].set_yscale('linear')

        ax[0].set_ylim([1, 2e4])
        ax[1].set_ylim([1e-4, 1e4])
        #ax[0,1].set_ylim([0.01, 1e5])
        #ax[1,0].set_ylim([1.01e-3, 1e2])

        #ax[0,0].set_xlim([1.01e-5, 0.99e2])
        #ax[0,1].set_xlim([1.01e-5, 0.99e2])
        #ax[1,0].set_xlim([1.01e-5, 0.99e2])
        #ax[1,1].set_xlim([1.01e-5, 0.99e2])

        ax[0].set_ylabel(r'Aerosol $\tau$ per bar', fontsize=24)
        ax[1].set_ylabel(r'Aerosol Cumulative $\tau$', fontsize=24)
        ax[2].set_ylabel('Temperature (K)', fontsize=24)

        ax[0].set_xlabel(r'Pressure (bar)', fontsize=24)
        ax[1].set_xlabel(r'Pressure (bar)', fontsize=24)
        ax[2].set_xlabel(r'Pressure (bar)', fontsize=24)

        #ax[2].axvline(x=ir_photosphere_pressure_bars,color='#c79fef', linestyle='dashed', linewidth=2)
        #ax[1].axvline(x=ir_photosphere_pressure_bars, color='#c79fef', linestyle='dashed',label='IR Photosphere', linewidth=2)
        #ax[0].axvline(x=ir_photosphere_pressure_bars, color='#c79fef', linestyle='dashed',  linewidth=2)

        ax[2].axvline(x=vis_photosphere_pressure_bars,color='#06c2ac', linestyle='dashed', linewidth=2)
        ax[1].axvline(x=vis_photosphere_pressure_bars, color='#06c2ac', linestyle='dashed',label='Optical Photosphere', linewidth=2)
        ax[0].axvline(x=vis_photosphere_pressure_bars, color='#06c2ac', linestyle='dashed',  linewidth=2)


        plt.rcParams['legend.title_fontsize'] = 'small'

        ax[1].yaxis.set_minor_locator(AutoMinorLocator(5))
        ax[1].tick_params(axis='y', which='minor', direction='in')

        ax[2].yaxis.set_minor_locator(AutoMinorLocator(5))
        ax[2].tick_params(axis='y', which='minor', direction='in')

        ax[0].legend(loc='best')
        ax[1].legend()
        ax[2].legend()

        fig.text(0.5, 0.92, r"Soot hazes, extended clouds, 100x metallicity", size=24, ha='center')

        plt.savefig('../Figures/Aerosol_Profiles_{}.png'.format(planet_name), bbox_inches='tight', dpi=100)
        i = i + 1
