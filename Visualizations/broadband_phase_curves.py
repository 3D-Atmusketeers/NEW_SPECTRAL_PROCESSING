#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as pl
import re
from scipy import interpolate

def record_phase_peak(phases, sw_fluxes, lw_fluxes, planet_name):
    """
    phases: list of phases or something
    fluxes: list of fluxes or something
    """

    num_points = 10000
    phases_new = np.linspace(0, 1, num_points)

    f_sw_fluxes = interpolate.interp1d(phases, sw_fluxes)
    f_lw_fluxes = interpolate.interp1d(phases, lw_fluxes)

    new_sw_fluxes = f_sw_fluxes(phases_new)
    new_lw_fluxes = f_lw_fluxes(phases_new)

    sw_peak = np.argmax(new_sw_fluxes) / num_points
    lw_peak = np.argmax(new_lw_fluxes) / num_points

    with open('OUTPUT_DATA/Phase_Curve_Peak_' + planet_name + '.txt', 'w') as f:
        f.write(planet_name + '\n')
        f.write("SW Peak: " + str(sw_peak) + "\n")
        f.write("LW PEAK: " + str(lw_peak) + "\n")
        f.write('\n\n')

    return None

def get_planet_df(planet_name, nlat, nlon, nlev, base):
    with open(base + planet_name + '/Planet_Run/fort.7', 'r') as fp:

        # read all lines using readline()
        lines = fp.readlines()
        for row in lines:

            # check if string present on a current line
            word = 'TOAALB'
            if row.find(word) != -1:
                TOAALB = float(re.findall(r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?", row)[0])

    with open(base + planet_name + '/Planet_Run/fort.7', 'r') as fp:

        # read all lines using readline()
        lines = fp.readlines()
        for row in lines:

            # check if string present on a current line
            word = 'FBASEFLUX'
            if row.find(word) != -1:
                FBASEFLUX = float(re.findall(r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?", row)[0])



    runfile = base + planet_name + '/Planet_Run/fort.62'

    # This is for the top of the weird column thing
    skip = 9
    latlon = np.empty([nlat * nlon, 2])  # to hold the latitude and longitude values at each location
    swout = np.empty([nlat * nlon])    # to hold the upward SW flux at each location (out the top of the atmosphere)
    swin = np.empty([nlat * nlon])    # to hold the upward SW flux at each location (out the top of the atmosphere)

    llind = np.linspace(nlev + 13, nlev + 13 + (nlat * nlon - 2) * (nlev + 11), num=nlat * nlon - 1)
    llextra = np.array([2])
    llind = np.concatenate((llextra, llind))                  # ugh, have to grab the right lines for lat, lon
    swind = np.linspace(skip, skip + (nlat * nlon - 1) * (nlev + 11), num=nlat * nlon)  # ugh, it's different for some reason here

    l = 0
    ind = 0
    with open(runfile) as f:
        for line in f:
            if l in llind:
                toss1, toss2, lat, lon = line.split()
                latlon[ind, :] = [np.float32(lat), np.float32(lon)]
            elif l in swind:
                sw, sw_in_val = line.split()

                rayleigh_reflected = TOAALB * float(sw_in_val)
                swout[ind] = float(sw) + rayleigh_reflected
                swin[ind] = float(sw_in_val) + rayleigh_reflected
                ind += 1
            l += 1
    f.close()


    # This is for the top of the weird column thing
    skip = 7
    lwout = np.empty([nlat * nlon])    # to hold the upward SW flux at each location (out the top of the atmosphere)

    llind   = np.linspace(nlev + 13, nlev + 13 + (nlat * nlon - 2) * (nlev + 11), num=nlat * nlon - 1)
    llextra = np.array([2])
    llind   = np.concatenate((llextra, llind))                  # ugh, have to grab the right lines for lat, lon
    swind   = np.linspace(skip, skip + (nlat * nlon - 1) * (nlev + 11), num=nlat * nlon)  # ugh, it's different for some reason here

    l = 0
    ind = 0
    with open(runfile) as f:
        for line in f:
            if l in llind:
                toss1, toss2, lat, lon = line.split()
                latlon[ind, :] = [np.float32(lat), np.float32(lon)]
            elif l in swind:
                lwout_val = float(line.split()[-1])
                lwout[ind] = lwout_val
                ind += 1
            l += 1
    f.close()

    lats = latlon[:, 0]
    lons = latlon[:, 1]

    df = pd.DataFrame({'lat': lats,
                       'lon': lons,
                       'swout': swout,
                       'swin': swin,
                       'lwout': lwout})

    df['lwin'] = FBASEFLUX

    df.lon = df.lon.mask(df.lon >= 180.0, df.lon - 360.0)

    return df


def get_phase_curve(df):
    df['swout_temp'] = df['swout'] * np.cos((df['lat'] * (np.pi / 180.0))) ** 2.0
    df['lwout_temp'] = df['lwout'] * np.cos((df['lat'] * (np.pi / 180.0))) ** 2.0
    df['swin_temp'] = df['swin'] * np.cos((df['lat'] * (np.pi / 180.0))) ** 2.0

    df['vis_temp'] = np.cos((df['lat'] * (np.pi / 180.0))) ** 2.0

    lon_df = df.groupby("lon").sum().reset_index()
    lon_df = lon_df.drop(['lat'], axis=1)

    # The arrays with the completed phase curves
    lw_phase_curve = []
    sw_phase_curve = []
    vis_phase_curve = []
    sw_in_phase_curve = []

    # For every lon, sum all the visible fluxes
    for i, lon in enumerate(list(np.unique(lon_df.lon))[::-1]):
        # Get the fraction that's normal to the observer for each lon
        lon_df['normal'] = np.cos((lon_df.lon - lon) * (np.pi / 180.))

        # Multiply the sw flux and the lw flux by that normal value
        lon_df['sw_final'] = lon_df['swout_temp'] * lon_df['normal']
        lon_df['lw_final'] = lon_df['lwout_temp'] * lon_df['normal']
        lon_df['vis_final'] = lon_df['vis_temp'] * lon_df['normal']

        lon_df['swin_final'] = lon_df['swin_temp'] * lon_df['normal']

        # Using the neat cos function, any value that's negative is not visible
        # Sum all the visible points
        lw_phase_curve.append(sum(filter(lambda x: x >= 0, list(lon_df.lw_final))))
        sw_phase_curve.append(sum(filter(lambda x: x >= 0, list(lon_df.sw_final))))
        sw_in_phase_curve.append(sum(filter(lambda x: x >= 0, list(lon_df.swin_final))))
        vis_phase_curve.append(sum(filter(lambda x: x >= 0, list(lon_df.vis_final))))

    # Divide by the visible area thing

    sw_curve_final = []
    lw_curve_final = []
    swin_curve_final = []
    for i in range(len(lw_phase_curve)):
        sw_curve_final.append(sw_phase_curve[i] / vis_phase_curve[i])
        lw_curve_final.append(lw_phase_curve[i] / vis_phase_curve[i])

        swin_curve_final.append(sw_in_phase_curve[i] / vis_phase_curve[i])

    return lw_curve_final, sw_curve_final, swin_curve_final, lon_df


def print_energy_balances(df, planet_name, sw_phase_curve, lw_phase_curve, nlon):
    df.swin = df.swin * np.cos((df['lat'] * (np.pi / 180.0)))
    df.lwin = df.lwin * np.cos((df['lat'] * (np.pi / 180.0)))

    df.swout = df.swout * np.cos((df['lat'] * (np.pi / 180.0)))
    df.lwout = df.lwout * np.cos((df['lat'] * (np.pi / 180.0)))

    total_lwout = sum(df.lwout)
    total_swout = sum(df.swout)
    total_lwin = sum(df.lwin)
    total_swin = sum(df.swin)

    orbital_phases = np.linspace(0, 1.0, nlon - 1)
    new_orbital_phases = np.linspace(0, 1.0, 100000)

    # Have to do :-1 because there is an extra 0
    sw_f = interpolate.interp1d(orbital_phases, sw_phase_curve[:-1], kind='quadratic')
    sw_interpolated_phase_curve = sw_f(new_orbital_phases)

    lw_f = interpolate.interp1d(orbital_phases, lw_phase_curve[:-1], kind='quadratic')
    lw_interpolated_phase_curve = lw_f(new_orbital_phases)

    with open('OUTPUT_DATA/Energy_Balances_' + planet_name + '.txt', 'w') as f:
        f.write('******************** \n')
        f.write(planet_name + '\n')
        f.write("Energy Out / Energy In: " + str(np.round((total_lwout + total_swout) / (total_lwin + total_swin), 3)) + "\n")
        f.write("BA: " + str(np.round(total_swout / total_swin, 3)) + "\n")
        f.write("Reflected Light Phase Curve Peak at: " + str(np.round(np.argmax(sw_interpolated_phase_curve) / 100000, 4)) + "\n")
        f.write("Thermal Light Phase Curve Peak at: " + str(np.round(np.argmax(lw_interpolated_phase_curve) / 100000, 4)) + "\n")
        f.write('\n\n')
    return None


def plot_thermal_phasecurves(planet_names, nlon, two_sets_of_planets):
    n = len(planet_names)
    colors = pl.cm.gnuplot(np.linspace(0, 1, int(n / 2) + 1))

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6), sharex=True, sharey=False)
    plt.subplots_adjust(wspace=0.275, hspace=0.1)

    for j, planet_name in enumerate(planet_names):
        z = j
        df = get_planet_df(planet_names[j], nlat=48, nlon=96, nlev=50, base='../Spectral-Processing/GCM-OUTPUT/')
        lw_phase_curve, sw_phase_curve, swin_phase_curve, lon_df = get_phase_curve(df)

        linestyle_str = 'solid'
        if two_sets_of_planets:
            if z >= int(n / 2):
                z = j - int(n / 2)
                linestyle_str = 'dashed'

        # Plot each phase curve
        temp_name = planet_names[j].replace("_", "-")
        phases = np.linspace(0, 1, nlon, endpoint=True)
        ax.plot(phases, lw_phase_curve, linewidth=2, label=temp_name, linestyle=linestyle_str, color=colors[z])

    ax.legend(fontsize=10, loc=(0, 1.05), ncol=2, mode='expand')

    ax.set_xlim(0.01, 0.99)
    ax.set_xlabel('Orbital Phase', fontsize=24)
    ax.set_xlabel('Orbital Phase', fontsize=24)
    ax.set_ylabel(r'Thermal Flux (W m$^{-2}$)', fontsize=24)
    plt.savefig('../Figures/broadband_thermal_phasecurves.png', bbox_inches='tight', dpi=100)
    plt.clf()

    return None


def plot_reflected_phasecurves(planet_names, nlon, two_sets_of_planets):
    n = len(planet_names)
    colors = pl.cm.gnuplot(np.linspace(0, 1, int(n / 2) + 1))

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6), sharex=True, sharey=False)
    plt.subplots_adjust(wspace=0.275, hspace=0.1)

    for j, planet_name in enumerate(planet_names):
        z = j
        df = get_planet_df(planet_names[j], nlat=48, nlon=96, nlev=50, base='../Spectral-Processing/GCM-OUTPUT/')
        lw_phase_curve, sw_phase_curve, swin_phase_curve, lon_df = get_phase_curve(df)

        linestyle_str = 'solid'
        if two_sets_of_planets:
            if z >= int(n / 2):
                z = j - int(n / 2)
                linestyle_str = 'dashed'

        temp_name = planet_names[j].replace("_", "-")
        phases = np.linspace(0, 1, nlon, endpoint=True)
        ax.plot(phases, sw_phase_curve, linewidth=2, label=temp_name, linestyle=linestyle_str, color=colors[z])

        record_phase_peak(phases, sw_phase_curve, lw_phase_curve, planet_name)
        print_energy_balances(df, planet_names[j], sw_phase_curve, lw_phase_curve, nlon)

    ax.legend(fontsize=10, loc=(0, 1.05), ncol=2, mode='expand')

    ax.set_xlim(0.01, 0.99)
    ax.set_xlim(0.01, 0.99)
    ax.set_xlabel('Orbital Phase', fontsize=24)
    ax.set_ylabel(r'Reflected Flux (W m$^{-2}$)', fontsize=24)

    plt.savefig('../Figures/broadband_reflected_phasecurves.png', bbox_inches='tight', dpi=100)
    plt.clf()

    return None


def plot_reflected_starlight_maps(planet_names):
    for j, planet_name in enumerate(planet_names):
        df = get_planet_df(planet_names[j], nlat=48, nlon=96, nlev=50, base='../Spectral-Processing/GCM-OUTPUT/')

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 5), sharex=True, sharey=False)
        plt.subplots_adjust(wspace=0.2, hspace=0.1)

        reflected_map = ax.scatter(df.lon, df.lat, c=df.swout, marker='s', s=30)
        ax.set_title(r'Reflected Flux (W m$^{-2}$)', fontsize=24)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.1)
        cbar1 = fig.colorbar(reflected_map, cax=cax, orientation='vertical')

        ax.set_xlim(-176, 176)
        ax.set_ylim(-87, 87)

        ax.set_xlabel('Longitude (degrees)', fontsize=24)
        ax.set_ylabel(r'Latitude (degrees)', fontsize=24)

        plt.savefig('../Figures/{}_reflected_starlight_map.png'.format(planet_names[j]), bbox_inches='tight', dpi=100)
        plt.clf()
    return None
