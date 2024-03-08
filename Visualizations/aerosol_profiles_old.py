import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import grab_input_data
import matplotlib.pylab as pl

def plot_aersol_profiles(planet_names, nlat, nlon, nlev, num_orders_of_magnitude):

    colors=['#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#ffffff']
    colors = pl.cm.nipy_spectral(np.linspace(0, 1, 13))

    for planet_name in planet_names:
        gravity = grab_input_data.get_input_data('../Spectral-Processing/GCM-OUTPUT/', planet_name,'fort.7' ,'GA')
        ir_absorbtion_coefficient = grab_input_data.get_input_data('../Spectral-Processing/GCM-OUTPUT/', planet_name,'fort.7' ,'ABSLW')

        ir_photosphere_pressure_bars = (2./3.) * (gravity/ir_absorbtion_coefficient) / 10000
        ir_photosphere_pressure_bars = np.round(ir_photosphere_pressure_bars, 3)


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
        
        fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(14,11))
        plt.subplots_adjust(wspace=0.25, hspace=0.15)

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


        if (np.max(cum_optical_depths_1) < 1e-10) and (np.max(cum_haze_optical_depths_1) < 1e-10):
            pass
        else:
            if 'ALL'.lower() in planet_name.lower():
                ax[0,0].plot(df1.pres, df1.aero_tau_2 / layer_pressures, color=colors[8], linewidth=3, label='ZnS')
                ax[0,0].plot(df1.pres, df1.aero_tau_3 / layer_pressures, color=colors[9], linewidth=3, label=r'Na$_2$S')
                ax[0,0].plot(df1.pres, df1.aero_tau_4 / layer_pressures, color=colors[10], linewidth=3, label='MnS')
                ax[0,0].plot(df1.pres, df1.aero_tau_9 / layer_pressures, color=colors[11], linewidth=3, label='Ni')
                ax[0,0].plot(df1.pres, df1.aero_tau_10 / layer_pressures, color=colors[12], linewidth=3, label='Fe')

                ax[0,1].plot(df2.pres, df1.aero_tau_2 / layer_pressures, color=colors[8], linewidth=3, label='ZnS')
                ax[0,1].plot(df2.pres, df1.aero_tau_3 / layer_pressures, color=colors[9], linewidth=3, label=r'Na$_2$S')
                ax[0,1].plot(df2.pres, df1.aero_tau_4 / layer_pressures, color=colors[10], linewidth=3, label='MnS')
                ax[0,1].plot(df2.pres, df1.aero_tau_9 / layer_pressures, color=colors[11], linewidth=3, label='Ni')
                ax[0,1].plot(df2.pres, df1.aero_tau_10 / layer_pressures, color=colors[12], linewidth=3, label='Fe')

            ax[0,0].plot(df1.pres, df1.aero_tau_1 / layer_pressures, color=colors[0], linewidth=3, label='KCl')
            #ax[0,0].plot(df1.pres, df1.aero_tau_5 / layer_pressures, color=colors[1], linewidth=3, label='Cr')
            #ax[0,0].plot(df1.pres, df1.aero_tau_6 / layer_pressures, color=colors[2], linewidth=3, label='SiO$_2$')
            #ax[0,0].plot(df1.pres, df1.aero_tau_7 / layer_pressures, color=colors[3], linewidth=3, label=r'Mg$_2$SiO$_4$')
            #ax[0,0].plot(df1.pres, df1.aero_tau_8 / layer_pressures, color=colors[4], linewidth=3, label='VO')
            #ax[0,0].plot(df1.pres, df1.aero_tau_11 / layer_pressures, color=colors[5], linewidth=3, label=r'Ca$_2$SiO$_4$')
            #ax[0,0].plot(df1.pres, df1.aero_tau_12 / layer_pressures, color=colors[6], linewidth=3, label=r'CaTiO3')
            #ax[0,0].plot(df1.pres, df1.aero_tau_13 / layer_pressures, color=colors[7], linewidth=3, label=r'Al$_2$O$_3$')
            ax[0,0].plot(df1.pres, haze_optical_depths_1 / layer_pressures, color='black', linewidth=3, linestyle='dashed', label=r'Haze')


            ax[0,1].plot(df2.pres, df2.aero_tau_1 / layer_pressures, color=colors[0], linewidth=3, label='KCl')
            #ax[0,1].plot(df2.pres, df2.aero_tau_5 / layer_pressures, color=colors[1], linewidth=3)
            #ax[0,1].plot(df2.pres, df2.aero_tau_6 / layer_pressures, color=colors[2], linewidth=3)
            #ax[0,1].plot(df2.pres, df2.aero_tau_7 / layer_pressures, color=colors[3], linewidth=3)
            #ax[0,1].plot(df2.pres, df2.aero_tau_8 / layer_pressures, color=colors[4], linewidth=3)
            #ax[0,1].plot(df2.pres, df2.aero_tau_11 / layer_pressures, color=colors[5], linewidth=3)
            #ax[0,1].plot(df2.pres, df2.aero_tau_12 / layer_pressures, color=colors[6], linewidth=3)
            #ax[0,1].plot(df2.pres, df2.aero_tau_13 / layer_pressures, color=colors[7], linewidth=3)
            ax[0,1].plot(df2.pres, haze_optical_depths_2 / layer_pressures, color='black', linewidth=3, linestyle='dashed',label=r'Haze')


            ax[1,0].plot(df1.pres, cum_optical_depths_1, color='black', linewidth=3, label='Clouds, Substellar Point')
            ax[1,0].plot(df2.pres, cum_optical_depths_2, color='red', linewidth=3, label='Clouds, Antistellar Point')

            ax[1,0].plot(df2.pres, cum_haze_optical_depths_1, color='black', linewidth=3,linestyle='dashed', label='Haze, Substellar Point')
            ax[1,0].plot(df2.pres, cum_haze_optical_depths_2, color='red', linewidth=3,linestyle='dashed',  label='Haze, Antistellar Point')


            ax[1,1].plot(df1.pres, df1.temp, color='black', linewidth=3, label='Substellar Profile')
            ax[1,1].plot(df2.pres, df2.temp, color='red', linewidth=3, label='Antistellar Profile')

            ax[0,0].set_xscale('log')
            ax[0,1].set_xscale('log')
            ax[1,0].set_xscale('log')
            ax[1,1].set_xscale('log')


            ax[0,0].set_yscale('log')
            ax[0,1].set_yscale('log')
            ax[1,0].set_yscale('log')
            ax[1,1].set_yscale('linear')

            ax[0,0].set_ylim([0.01, 1e5])
            ax[0,1].set_ylim([0.01, 1e5])
            ax[1,0].set_ylim([1.01e-3, 1e2])

            ax[0,0].set_xlim([1.01e-5, 0.99e2])
            ax[0,1].set_xlim([1.01e-5, 0.99e2])
            ax[1,0].set_xlim([1.01e-5, 0.99e2])
            ax[1,1].set_xlim([1.01e-5, 0.99e2])

            ax[0,0].set_ylabel(r'Aerosol $\tau$ per bar', fontsize=20)
            ax[0,1].set_ylabel(r'Aerosol $\tau$ per bar', fontsize=20)
            ax[1,0].set_ylabel(r'Aerosol Cumulative $\tau$', fontsize=20)
            ax[1,1].set_ylabel('Temperature (K)', fontsize=20)

            ax[0,0].axvline(x=ir_photosphere_pressure_bars,color='gray', linestyle='dashed', linewidth=3)
            ax[0,1].axvline(x=ir_photosphere_pressure_bars, color='gray', linestyle='dashed', linewidth=3)
            ax[1,0].axvline(x=ir_photosphere_pressure_bars, color='gray', linestyle='dashed', linewidth=3)
            ax[1,1].axvline(x=ir_photosphere_pressure_bars, color='gray', linestyle='dashed', label='Photosphere Pressure Level', linewidth=3)

            plt.rcParams['legend.title_fontsize'] = 'small'

                        
            #ax.legend(fontsize=12, loc=(0, 1.05), ncol=2, mode='expand', title_fontsize=16)
            
                        
            #plt.figlegend(loc=(0, 1.05), ncol=7, labelspacing=0., title="Cloud Species", fontsize=14)

                        
            #ax[0,0].legend(fontsize=12, loc=(0, 1.05), ncol=6, mode='expand', title="Cloud Species", title_fontsize=16)            
            ax[1,0].legend(fontsize=12)
            ax[1,1].legend(fontsize=12)
            ax[0,0].legend(fontsize=12, ncol=1, labelspacing=0.0, loc='upper left')
            ax[0,1].legend(fontsize=12, ncol=1, labelspacing=0.0, loc='upper left')

            fig.text(0.39, 0.85, r"Substellar Profile", size=14, ha='center')
            fig.text(0.82, 0.85, r"Antistellar Profile", size=14, ha='center')

            fig.text(0.5, 0.04, r"Pressure (bar)", size=20, ha='center')

            plt.savefig('../Figures/Aerosol_Profiles_{}.png'.format(planet_name), bbox_inches='tight', dpi=100)
            i = i + 1