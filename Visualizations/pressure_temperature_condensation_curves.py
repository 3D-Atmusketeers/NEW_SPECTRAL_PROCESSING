#!/usr/bin/env python
# coding: utf-8

### ----- IMPORT LIBRARIES ----- ###
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams, rc
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
import math

def plot_PTC_curves(planet_names, nlat, nlon, nlev, num_gcms, nucleation_lim):
    base = '../Spectral-Processing/PLANET_MODELS/'

    column_names = ['lat', 'lon', 'level',
                'alt', 'pres', 'temp', 
                'u', 'v', 'w']

    nparams = len(column_names)

    ### ----- INPUT/OUTPUT CONTROL ----- ###

    condensation_data = 'DATA/Condensation_Ts.dat'       # condensation curves

    Tmin = 0
    Tmax = 650

    NPARAMS = 9

    for ind, planet_name in enumerate(planet_names):
        plt.clf()
        print('Creating plot for', planet_name)
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
        plt.subplots_adjust(wspace=0.02, hspace=0.03)

        if 'CLOUDS' in planet_name:
            clouds = True
        else:
            clouds = False

        df = pd.read_csv(base + planet_name + '.txt', names = column_names, delim_whitespace=True)
        df.lon = df.lon.mask(df.lon >= 180.0, df.lon - 360.0)
        max_temp = math.ceil(df['temp'].max() / 100) * 100
        min_temp = max(math.floor(df['temp'].min() / 100) * 100, 0)

        
        data = np.loadtxt(base + planet_name + '.txt')
        data = data.reshape(nlat, nlon, nlev, NPARAMS)

        lons = data[:,:,0,1][0]
        lats = data[:,:,0,0][:,0]
        
        delta_temps = np.zeros((nlat, nlon))
        taus = np.zeros((nlat, nlon))


        # colormap
        cm_name = 'nuuk'
        cm_file = np.loadtxt(f'ScientificColourMaps7/{cm_name}/{cm_name}.txt')
        my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file)

        # adjust tick marks on both axes

        axes.tick_params(axis='both',
                    which='both',
                    direction='in',
                    top = True,
                    right = True)
        axes.tick_params(axis='both',
                    which='major',
                    length=5)


        # plot all TP-profiles (expensive)
        for i in range(48):
            for j in range(96):
                axes.semilogy(data[i][j][:,5], data[i][j][:,4], alpha=0.2, color='gray', linewidth=1)


        # colormap
        #cm_name = 'batlow'
        cm_name = 'bamO'
        cm_file = np.loadtxt(f'ScientificColourMaps7/{cm_name}/{cm_name}.txt')
        cm_file  = np.roll(cm_file, 140, axis=0)
        my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file)

        colors = np.linspace(0, 1, 96)
        # set color cycle
        color_idx = 48
        
        # plot colored equatorial TP-profiles
        for i in range(96):
            axes.semilogy(data[24][i][:,5], data[24][i][:,4],
                                '-', lw=1.8, alpha=0.75, color=my_colors(colors[color_idx]))

            color_idx += 1
            if (color_idx > 95):
                color_idx = color_idx - 96
        
        # plot cloud condensation curves
        if clouds == True:
            press_con, kcl, zns, na2s, mns, cr, sio2, mg2sio4, vo, ni, fe, ca2sio4, catio3, al2o3, mgsio3= np.loadtxt(condensation_data, unpack=True)
            # KCl, Cr, SiO2, Mg$_2$SiO$_4$, VO, Ni, Ca$_2$SiO$_4$, CaTiO$_3$, Al$_2$O$_3$
            
            new_pressures = np.logspace(-5, 2, 200)
            
            f_kcl = interp1d( press_con, kcl, fill_value='extrapolate')
            f_zns = interp1d(press_con, zns,  fill_value='extrapolate')
            f_na2s = interp1d(press_con,na2s,  fill_value='extrapolate')
            f_mns = interp1d(press_con, mns, fill_value='extrapolate')
            f_cr = interp1d(press_con,cr,  fill_value='extrapolate')
            f_sio2 = interp1d(press_con,sio2,  fill_value='extrapolate')
            f_mg2sio4 = interp1d(press_con,mg2sio4,  fill_value='extrapolate')
            f_vo = interp1d(press_con,vo,  fill_value='extrapolate')
            f_ni = interp1d(press_con,ni,  fill_value='extrapolate')
            f_fe = interp1d(press_con,fe,  fill_value='extrapolate')
            f_ca2sio4 = interp1d(press_con, ca2sio4,  fill_value='extrapolate')
            f_catio3 = interp1d(press_con,catio3,  fill_value='extrapolate')
            f_al2o3 = interp1d( press_con,al2o3, fill_value='extrapolate')
            
            colors=['#4b006e', '#3f9b0b', '#75bbfd', '#ff81c0', 'brown', '#fb2943', '#029386', 'orange']

            axes.semilogy(f_kcl(new_pressures), new_pressures, lw=3, label='KCl', color=colors[0])
            axes.semilogy(f_cr(new_pressures), new_pressures, lw=3, label='Cr', color=colors[1])
            axes.semilogy(f_sio2(new_pressures), new_pressures, lw=3, label=r'SiO$_2$', color=colors[2])
            axes.semilogy(f_mg2sio4(new_pressures), new_pressures, lw=3, label=r'Mg$_2$SiO$_4$', color=colors[3])
            axes.semilogy(f_vo(new_pressures), new_pressures, lw=3, label='VO', color=colors[4])
            axes.semilogy(f_ca2sio4(new_pressures), new_pressures, lw=3, label=r'Ca$_2$SiO$_4$', color=colors[5])
            axes.semilogy(f_catio3(new_pressures), new_pressures, lw=3, label=r'CaTiO$_3$', color=colors[6])
            axes.semilogy(f_al2o3(new_pressures), new_pressures, lw=3, label=r'Al$_2$O$_3$', color=colors[7])


            if nucleation_lim == False:
                axes.semilogy(f_fe(new_pressures), new_pressures, lw=3, label=r'Fe', color='darkgreen')
                axes.semilogy(f_na2s(new_pressures), new_pressures, lw=3, label=r'Na$_2$S', color='mediumseagreen')
                axes.semilogy(f_mns(new_pressures), new_pressures, lw=3, label=r'MnS', color='lightseagreen')
                axes.semilogy(f_zns(new_pressures), new_pressures, lw=3, label=r'ZnS', color='teal')
                axes.semilogy(f_ni(new_pressures), new_pressures, lw=3, label=r'Ni', color='steelblue')
                
            axes.legend(fontsize=14, ncol=2, labelspacing=0.0, loc='lower left')



        
        axes.set_ylim([1.01e-5, 1.00e+2])
        axes.set_xlim([min_temp,max_temp])

        
        axes.invert_yaxis()
        axes.xaxis.set_ticks_position('bottom')
        axes.xaxis.set_label_position('bottom')
        
        
        # tp labels
        axes.set_ylabel('Pressure (bars)', fontsize=24)


        # colorbars
        sm = plt.cm.ScalarMappable(cmap=my_colors, norm=plt.Normalize(vmin=-180, vmax=180))
        sm._A = []
        tp_cbar = fig.colorbar(sm, aspect=15, pad=0.02, ticks=[0, -60, -120, 60, 120])
        tp_cbar.set_label('Lon (deg)', fontsize=24, labelpad=5)


        fig.text(0.5, 0.95, r"Temperature (K)", size=24, ha='center')


        print('Creating plot DONE')

        fig.savefig('../Figures/PTC_Curves_{}.png'.format(planet_name), bbox_inches='tight', dpi=100)
