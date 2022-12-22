#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import numpy as np
import math
from scipy.interpolate import griddata
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as patches
import matplotlib as mpl
from matplotlib.font_manager import FontProperties
import numpy as np
import math
from scipy.io import readsav
from astropy import units as u
from astropy import constants as const
import matplotlib.pylab as pl
from scipy import interpolate
from scipy import ndimage
import numpy as np
from scipy import ndimage as ndi
import scipy.interpolate as spi
from scipy.integrate import trapz
from astropy.modeling import models
from astropy import units as u
from astropy.modeling.models import BlackBody
from astropy.visualization import quantity_support
import re


# GLOBAL PARAMS
h = 6.6260701e-34
c = 3e8

STAR_RADIUS   = 0.215 * 695700000
PLANET_RADIUS = 1.7469282e7

def get_filter(which_filter):
    #df_5 = {'wav': np.linspace(5e-6, 6e-6, 500), 'transmission': np.linspace(1, 1, 500)}
    

    #filter_3_6 = pd.read_csv('DATA/filter_3.6.dat', delim_whitespace=True, skiprows=3, names=['wav', 'transmission'])
    #filter_4_5 = pd.read_csv('DATA/filter_4.5.dat', delim_whitespace=True, skiprows=3, names=['wav', 'transmission'])
    #filter_5_0 = pd.DataFrame(data=df_5)
    

    # Put the filters in meters
    #filter_3_6.wav = filter_3_6.wav / 1e6
    #filter_4_5.wav = filter_4_5.wav / 1e6
    

    # Interpolate the filters
    #f_3_6 = interpolate.interp1d(filter_3_6.wav, filter_3_6.transmission)
    #f_4_5 = interpolate.interp1d(filter_4_5.wav, filter_4_5.transmission)
    #f_5_0 = interpolate.interp1d(filter_5_0.wav, filter_5_0.transmission)


    if which_filter == 'MIRI':
        df_5_12 = {'wav': np.linspace(5e-6, 12e-6, 500), 'transmission': np.linspace(1, 1, 500)}
        filter_5_12 = pd.read_csv('DATA/MIRI_BANDPASS.txt', delim_whitespace=True, skiprows=1, names=['wav', 'transmission'])
        filter_5_12.wav = filter_5_12.wav / 1e6
        f_5_12 = interpolate.interp1d(filter_5_12.wav, filter_5_12.transmission)
        return filter_5_12, f_5_12
    else:
        print ("Filters other than MIRI are not written in yet")
        exit(0)
        

    


def get_star_spectra():
    # Get the star spectra stuff
    stellar_spectrum_1 = pd.read_csv('DATA/model_3250K.txt', names=['wavelength', 'flux_Watt/m^2/micron'], delim_whitespace=True, skiprows=1)
    wavelengths_meters    = np.asarray(list(stellar_spectrum_1.wavelength * 1e-6)) * u.m
    flux_Watt_m_2_microns = np.asarray(list(stellar_spectrum_1['flux_Watt/m^2/micron'])) * u.W / u.m**2 / u.um
    flux_si = flux_Watt_m_2_microns * np.pi
    star_spectra = pd.DataFrame(list(zip(wavelengths_meters.value,flux_si.value)), columns=['wavelength','flux'])
    return star_spectra


def specint(wn, spec, filtwn_list, filttrans_list): #from Ryan Challener's Theresa code
    """
    Integrate a spectrum over the given filters.
    Arguments
    ---------
    wn: 1D array
        Wavenumbers (/cm) of the spectrum 
        Hayley's note: be wary of these units! 
    spec: 1D array
        Spectrum to be integrated
    filtwn_list: list
        List of arrays of filter wavenumbers, in /cm.
    filttrans_list: list
        List of arrays of filter transmission. Same length as filtwn_list.
        Hayley's note: for now, this can just be an array you create of 1's and 0's. (1 where want to integrate over, 0 outside)
    Returns
    -------
    intspec: 1D array
        The spectrum integrated over each filter. 
    """
    if len(filtwn_list) != len(filttrans_list):
        print("ERROR: list sizes do not match.")
        raise Exception
    
    intspec = np.zeros(len(filtwn_list)) 
    
    for i, (filtwn, filttrans) in enumerate(zip(filtwn_list, filttrans_list)):
        # Sort ascending
        idx = np.argsort(filtwn)
        
        intfunc = spi.interp1d(filtwn[idx], filttrans[idx],
                               bounds_error=False, fill_value=0)

        # Interpolate transmission
        inttrans = intfunc(wn)

        # Normalize to one
        norminttrans = inttrans / np.trapz(inttrans, wn)
        

        # Integrate filtered spectrum
        intspec[i] = np.trapz(spec * norminttrans, wn)

    return intspec


def getSpacing(arr):
    return (arr[-1]-arr[0])/(len(arr)-1)

def fwhm2sigma(fwhm):
    return fwhm / np.sqrt(8 * np.log(2))

def gaussian(x, mu, sig):
    normalization = 1/np.sqrt(2 * np.pi * sig**2)
    exponent = - ((x - mu)**2 / (2*sig**2))
    y = normalization * np.exp(exponent)
    return y

def reduceSpectralResolution(x, y, R_low, R_high=None, lambda_mid=None, n=4):
    if np.ndim(y) == 2:
        return np.array([reduceSpectralResolution(x,spec,R_low,R_high,lambda_mid,n) for spec in y])

    dx = getSpacing(x)
    # If lambda_mid is none, take median of input wavelengths
    if lambda_mid is None:
        lambda_mid = np.median(x)

    # If R_high is none, use midpoint of x divided by spacing in x
    if R_high is None:
        R_high = lambda_mid/dx

    # Create Gaussian kernel
    fwhm = np.sqrt(R_high**2 - R_low**2)*(lambda_mid/(R_low*R_high))
    sigma = fwhm2sigma(fwhm)

    kernel_x = np.arange(-n*sigma, n*sigma+dx, dx)
    kernel = gaussian(kernel_x,0,sigma)
    kernel = kernel/np.sum(kernel)

    # find center of kernel
    n_kernel_lt0 = len(np.where(kernel_x<0)[0])
    n_kernel_gt0 = len(np.where(kernel_x>0)[0])

    if n_kernel_lt0< n_kernel_gt0:
        origin = 0
    else:
        origin=-1

    # convolve
    lowRes = ndimage.convolve(y, kernel, origin=origin)
    return lowRes


# ## Spectra Test
def plot_planet_spectra_blackbody_comparison(planet_names, black_body_temperatures, num_phases):
    """
    plot the planet spectra, compared to a blackbody
    planet names: list of strings
    black_body_temperatures: list of floats
    """

    colors_bb = pl.cm.Greys(np.linspace(0,1,len(black_body_temperatures) + 2))

    cm_name = 'batlow'
    cm_file = np.loadtxt(f'ScientificColourMaps7/{cm_name}/{cm_name}.txt')
    cm_file  = np.roll(cm_file, 140, axis=0)
    my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file)


    for k, planet_name in enumerate(planet_names):
        # Figure aesthetics
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6),sharex=True, sharey=False)
        plt.subplots_adjust(hspace=0.05, wspace=0.25)


        for i in range(num_phases):
            # Get the phase value
            rot_val = 360. / num_phases
            
            # Get the file path
            file_path = '../Spectral-Processing/FINISHED_SPECTRA/Spec_0_' + planet_name + '_phase_{}_inc_0.00.00.0000.00.dat'
            
            # Load in the planet spectra
            planet_spectra = pd.read_csv(file_path.format(str(i * rot_val)), header=None, delim_whitespace=True, names=['wavelength','flux', 'reflected'])
    
            # Reset the index
            planet_spectra = planet_spectra.reset_index(drop=True)

            ax.plot(planet_spectra.wavelength * 1e6,
                    planet_spectra.flux,
                    color=my_colors(i / num_phases),
                    linewidth=3,
                    label="Planet Spectra at phase: " + str(np.round(rot_val * i / 360., 3)))

            # Plot the blackbody temperatures
            j = 0
            for temp in black_body_temperatures:
                bb = BlackBody(temperature=temp*u.K)
                wav = np.linspace(0.5, 20.0, 1000) * u.um
                flux = bb(wav) * np.pi * u.sr

                # Convert to si, same as the output of the emission code
                flux = flux.to(u.J / (u.m * u.m * u.Hz * u.s))

                ax.plot(wav, flux,
                        color=colors_bb[j + 2],
                        alpha=0.8,
                        linewidth=2, 
                        linestyle='dashed',
                        label='Blackbody at: ' + str(temp) + 'K')
                j = j + 1
        
        # Center the y value
        center_y_val = int(np.log10(np.mean(planet_spectra.flux)))

        # Do somet figure stuff
        ax.set_ylim(10 ** (center_y_val - 2), 10 ** (center_y_val + 2))
        ax.set_yscale('log')
        ax.set_xlim(min(planet_spectra.wavelength*1e6),max(planet_spectra.wavelength*1e6))          
        ax.legend(fontsize=12, loc=(0,1.05), ncol=2, mode='expand', title_fontsize=16)
        ax.set_xlabel('Wavelength ($\mu$m)')
        ax.set_ylabel(r'Flux (W/m$^2 Hz$)') #  (W m$^{-2}$)
        plt.savefig('../Figures/planet_spectra_blackbody_comparison_{}.jpg'.format(planet_name), dpi=200, bbox_inches='tight')
        plt.clf()        
    return None


# ## Plot the star

# In[ ]:

def plot_star_spectra_test(planet_names):
    # Figure aesthetics
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6),sharex=True, sharey=False)
    plt.subplots_adjust(hspace=0.05, wspace=0.25)

    # get the star_spectra
    star_spectra = get_star_spectra()  

    plt.plot(star_spectra.wavelength * 1e6, star_spectra.flux, label='Stellar Spectrum', alpha=0.8, color='black', linewidth=2)

    plt.ylim(1, 1e7)
    plt.yscale('log')
    plt.xlabel('Wavelengths (microns)')
    plt.ylabel(r'Flux (Watts/m$^2$/micron)')
    plt.xlim(0, 6)
    plt.legend(loc='upper right')
    plt.savefig('../Figures/star_test.jpg', dpi=200, bbox_inches='tight')
    plt.clf()

    return None


def plot_filters(planet_names):
    filter_5_12, f_5_12 = get_filter(which_filter='MIRI')

    # Figure aesthetics
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
    plt.subplots_adjust(hspace=0.05, wspace=0.25)

    # Load in a test spectra to get the wavelength ranges
    spectra = pd.read_csv('../Spectral-Processing/FINISHED_SPECTRA/Spec_0_GJ1214b-HAZES-100X-soot_phase_0.0_inc_0.00.00.0000.00.dat', 
                        header=None, delim_whitespace=True, names=['wavelength','flux', 'reflected'])

    # Get only the spectra in the filter ranges we care about
    #spectra_3_6  = spectra[(spectra['wavelength'] > np.amin(filter_3_6.wav)) & (spectra['wavelength'] < np.amax(filter_3_6.wav))]  
    #spectra_4_5  = spectra[(spectra['wavelength'] > np.amin(filter_4_5.wav)) & (spectra['wavelength'] < np.amax(filter_4_5.wav))]  
    spectra_5_12 = spectra[(spectra['wavelength'] > np.amin(filter_5_12.wav)) & (spectra['wavelength'] < np.amax(filter_5_12.wav))]  


    # Plot the interpolated transmission stuff and the original ones
    #ax.plot(filter_3_6.wav, filter_3_6.transmission, color='#7e1e9c', linewidth=3, alpha=0.2, label="")
    #ax.plot(filter_4_5.wav, filter_4_5.transmission, color='#15b01a', linewidth=3, alpha=0.2, label="Original")
    ax.plot(filter_5_12.wav * 1e6, filter_5_12.transmission, color='black', linewidth=3, alpha=0.2, label="MIRI")

    #ax.plot(spectra_3_6.wavelength, f_3_6(spectra_3_6.wavelength),color='#7e1e9c', linewidth=3, linestyle='dashed', label="Interpolated")
    #ax.plot(spectra_4_5.wavelength, f_4_5(spectra_4_5.wavelength),color='#15b01a', linewidth=3, linestyle='dashed', label="Interpolated")
    ax.plot(spectra_5_12.wavelength * 1e6, f_5_12(spectra_5_12.wavelength),color='black', linewidth=3, linestyle='dashed', label="Interpolated")

    # Some more Figure stuff
    ax.set_xlabel(r"Wavelength ($\mu$m)")
    ax.set_ylabel('Transmittance')
    ax.legend(loc='lower right', ncol=2, fontsize=14)
    plt.savefig('../Figures/Filters.jpg', dpi=100, bbox_inches='tight')
    plt.clf()
    return None

def plot_spectra_simple(planet_names, num_phases):
    cm_name = 'batlow'
    cm_file = np.loadtxt(f'ScientificColourMaps7/{cm_name}/{cm_name}.txt')
    cm_file  = np.roll(cm_file, 140, axis=0)
    my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file)

    for k, planet_name in enumerate(planet_names):
        # Figure aesthetics
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6),sharex=True, sharey=False)
        plt.subplots_adjust(hspace=0.05, wspace=0.25)


        for i in range(num_phases):
            # Get the phase value
            rot_val = 360. / num_phases
            
            # Get the file path
            file_path = '../Spectral-Processing/FINISHED_SPECTRA/Spec_0_' + planet_name + '_phase_{}_inc_0.00.00.0000.00.dat'
            
            # Load in the planet spectra
            planet_spectra = pd.read_csv(file_path.format(str(i * rot_val)), header=None, delim_whitespace=True, names=['wavelength','flux', 'reflected'])
    
            # Reset the index
            planet_spectra = planet_spectra.reset_index(drop=True)

            ax.plot(planet_spectra.wavelength * 1e6,
                    planet_spectra.flux,
                    color=my_colors(i / num_phases),
                    linewidth=3,
                    label=str(np.round(rot_val * i / 360., 3)))
        
        # Do somet figure stuff
        ax.set_xlim(min(planet_spectra.wavelength*1e6),max(planet_spectra.wavelength*1e6))          
        ax.legend(fontsize=12, loc=(0,1.08), ncol=5, mode='expand', title='Orbital Phase', title_fontsize=16)
        ax.set_xlabel('Wavelength ($\mu$m)')
        ax.set_ylabel(r'Flux (W/m$^2$/$\mu m$)') #  (W m$^{-2}$)
        plt.savefig('../Figures/Planet_Simple_Spectra_{}.jpg'.format(planet_name), dpi=200, bbox_inches='tight')
        plt.clf()

    return None


def plot_spectra_phases(planet_names, num_phases, transmission_filter_name, planet_only_bool):
    cm_name = 'batlow'
    cm_file = np.loadtxt(f'ScientificColourMaps7/{cm_name}/{cm_name}.txt')
    cm_file  = np.roll(cm_file, 140, axis=0)
    my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file)
     
    # I currently only have MIRI coded in
    if transmission_filter_name == 'MIRI':
        filt, interp_function = get_filter(which_filter=transmission_filter_name)
    else:
        pass

    for k, planet_name in enumerate(planet_names):
        # Figure aesthetics
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6),sharex=True, sharey=False)
        plt.subplots_adjust(hspace=0.05, wspace=0.25)

        integrated_signal     = []
        integrated_signal_star = []

        for i in range(num_phases):
            rot_val = 360. / num_phases
            integrated_signal.append(0)
            integrated_signal_star.append(0)

            file_path = '../Spectral-Processing/FINISHED_SPECTRA/Spec_0_' + planet_name + '_phase_{}_inc_0.00.00.0000.00.dat'
            
            # Load in the planet spectra and convert to per microns
            # The star spectrum is per micron
            planet_spectra = pd.read_csv(file_path.format(str(i * rot_val)), header=None, delim_whitespace=True, names=['wavelength','flux', 'reflected'])
            planet_spectra.flux = planet_spectra.flux * (3.0e8 / planet_spectra.wavelength ** 2) / 1e6
    
            # Get the star spectra
            star_spectra = get_star_spectra()   
            
            if transmission_filter_name == 'MIRI':
                # Only use the data points within the filter wavelength range that we want
                planet_spectra  = planet_spectra[(planet_spectra['wavelength'] > np.amin(filt.wav)) & (planet_spectra['wavelength'] < np.amax(filt.wav))] 
                star_spectra = star_spectra[(star_spectra['wavelength'] > np.amin(filt.wav)) & (star_spectra['wavelength'] < np.amax(filt.wav))]  
            else:
                pass

            # Reset the indexes of the planet and star dataframes
            planet_spectra = planet_spectra.reset_index(drop=True)
            #star_spectra   = star_spectra.reset_index(drop=True)

            if planet_only_bool:
                # Convolve the planet and star down to a lower resolution
                reduced_planet = np.asarray(reduceSpectralResolution(list(planet_spectra.wavelength),list(planet_spectra.flux),R_low=100, R_high=None, lambda_mid=None, n=4))

                ax.plot(planet_spectra.wavelength * 1e6,
                        reduced_planet,
                        color=my_colors(i / num_phases),
                        linewidth=3,
                        label=str(np.round(rot_val * i / 360., 3)))

                pd.DataFrame({'Wave (microns)':planet_spectra.wavelength * 1e6,
                        'Planet_Flux SI/micron':reduced_planet}).to_csv('OUTPUT_DATA/Planet_White_Light_Phase_{}_Spectra_{}.txt'.format(str(i * rot_val), planet_name), sep=' ')


            else:
                # Convolve the planet and star down to a lower resolution
                reduced_planet = np.asarray(reduceSpectralResolution(list(planet_spectra.wavelength),
                                                                    list(planet_spectra.flux),
                                                                    R_low=100, R_high=None, lambda_mid=None, n=4))
                reduced_star = np.asarray(reduceSpectralResolution(list(star_spectra.wavelength),
                                                                list(star_spectra.flux),
                                                                R_low=100, R_high=None, lambda_mid=None, n=4))

                # Regrid the star to the same wavelength grid as the planet
                fnc_star_converter      = interpolate.interp1d(star_spectra.wavelength, reduced_star, fill_value='extrapolate')
                regridded_star_spectrum = fnc_star_converter(planet_spectra.wavelength)


                # Calculate the total signal, and the star signal
                if transmission_filter_name == 'MIRI':
                    signal      = (regridded_star_spectrum * STAR_RADIUS**2. + reduced_planet * PLANET_RADIUS**2.) * (planet_spectra.wavelength / h*c) * interp_function(planet_spectra.wavelength)
                    signal_star = (regridded_star_spectrum * STAR_RADIUS**2.) * (planet_spectra.wavelength / h*c) * interp_function(planet_spectra.wavelength)
                else:
                    signal      = (regridded_star_spectrum * STAR_RADIUS**2. + reduced_planet * PLANET_RADIUS**2.) * (planet_spectra.wavelength / h*c) * 1.
                    signal_star = (regridded_star_spectrum * STAR_RADIUS**2.) * (planet_spectra.wavelength / h*c)  * 1.


                # Divide out the two integrated signals
                fp_fs_ratio = (signal / signal_star - 1.0)

                ax.plot(planet_spectra.wavelength * 1e6, fp_fs_ratio * 1e6,
                        color=my_colors(i / num_phases),
                        linewidth=3,
                        label=str(np.round(rot_val * i / 360., 3)))
            
    ax.set_xlim(min(planet_spectra.wavelength*1e6),max(planet_spectra.wavelength*1e6))  
    ax.legend(fontsize=12, loc=(0,1.08), ncol=5, mode='expand', title='Orbital Phase', title_fontsize=18)
    ax.set_xlabel('Wavelength ($\mu$m)')

    if planet_only_bool:
        ax.set_ylabel(r'Planet Flux (W/m$^2$/$\mu m$)') #  (W m$^{-2}$)
        plt.savefig('../Figures/Planet_White_Light_Spectra_{}.jpg'.format(planet_name), dpi=200, bbox_inches='tight')
    else:
        ax.set_ylabel(r'F$_p$/F$_s$ (ppm)') #  (W m$^{-2}$)
        plt.savefig('../Figures/Ratio_White_Light_Spectra_{}.jpg'.format(planet_name), dpi=200, bbox_inches='tight')
    return None


def plot_phase_curves(planet_names, planet_name_char_len,  num_phases, transmission_filter_name):
    # Figure aesthetics
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6),sharex=True, sharey=False)
    plt.subplots_adjust(hspace=0.05, wspace=0.25)


    cm_name = 'batlow'
    cm_file = np.loadtxt(f'ScientificColourMaps7/{cm_name}/{cm_name}.txt')
    cm_file  = np.roll(cm_file, 140, axis=0)
    my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file)
    
    # I currently only have MIRI coded in
    if transmission_filter_name == 'MIRI':
        filt, interp_function = get_filter(which_filter=transmission_filter_name)
    else:
        pass

    for k, planet_name in enumerate(planet_names):
        integrated_signal     = []
        integrated_signal_star = []

        
        for i in range(num_phases):
            rot_val = 360. / num_phases
            integrated_signal.append(0)
            integrated_signal_star.append(0)
            
            # Load in the planet spectra and convert to per microns
            # The star spectrum is per micron
            
            file_path = '../Spectral-Processing/FINISHED_SPECTRA/Spec_0_' + planet_name + '_phase_{}_inc_0.00.00.0000.00.dat'
            planet_spectra = pd.read_csv(file_path.format(str(i * rot_val)), header=None, delim_whitespace=True, names=['wavelength','flux', 'reflected'])
            planet_spectra.flux = planet_spectra.flux * (3.0e8 / planet_spectra.wavelength ** 2) / 1e6
                    
            # Get the star spectra
            star_spectra = get_star_spectra()   
  
            if transmission_filter_name == 'MIRI':
                planet_spectra  = planet_spectra[(planet_spectra['wavelength'] > np.amin(filt.wav)) & (planet_spectra['wavelength'] < np.amax(filt.wav))] 
                star_spectra    = star_spectra[(star_spectra['wavelength'] > np.amin(filt.wav)) & (star_spectra['wavelength'] < np.amax(filt.wav))]  

                # Reset the indexes of the planet and star dataframes
                planet_spectra = planet_spectra.reset_index(drop=True)
                star_spectra   = star_spectra.reset_index(drop=True)

            else:
                pass

            # Regrid the star to the same wavelength grid as the planet
            fnc_star_converter      = interpolate.interp1d(star_spectra.wavelength, star_spectra.flux, fill_value='extrapolate')
            regridded_star_spectrum = fnc_star_converter(planet_spectra.wavelength)
            
            if transmission_filter_name == 'MIRI':
                signal      = (regridded_star_spectrum * STAR_RADIUS**2. + planet_spectra.flux * PLANET_RADIUS**2.) \
                            * (planet_spectra.wavelength / h*c) * interp_function(planet_spectra.wavelength)

                signal_star = (regridded_star_spectrum * STAR_RADIUS**2.) \
                            * (planet_spectra.wavelength / h*c) * interp_function(planet_spectra.wavelength)
            else:
                signal      = (regridded_star_spectrum * STAR_RADIUS**2. + planet_spectra.flux * PLANET_RADIUS**2.) \
                        * (planet_spectra.wavelength / h*c) * 1.

                signal_star = (regridded_star_spectrum * STAR_RADIUS**2.) \
                        * (planet_spectra.wavelength / h*c) * 1.
            
            # Convolve the planet and star down to a lower resolution
            #reduced_planet = np.asarray(reduceSpectralResolution(list(planet_spectra.wavelength),
            #                                                     list(planet_spectra.flux),
            #                                                     R_low=100, R_high=None, lambda_mid=None, n=4))
            #reduced_star = np.asarray(reduceSpectralResolution(list(star_spectra.wavelength),
            #                                                   list(star_spectra.flux),
            #                                                   R_low=100, R_high=None, lambda_mid=None, n=4))
            
            # Regrid the star to the same wavelength grid as the planet
            fnc_star_converter      = interpolate.interp1d(star_spectra.wavelength, star_spectra.flux, fill_value='extrapolate')
            regridded_star_spectrum = fnc_star_converter(planet_spectra.wavelength)
    
        
            if transmission_filter_name == 'MIRI':
                signal      = (regridded_star_spectrum * STAR_RADIUS**2. + planet_spectra.flux * PLANET_RADIUS**2.) \
                            * (planet_spectra.wavelength / h*c) * interp_function(planet_spectra.wavelength)

                signal_star = (regridded_star_spectrum * STAR_RADIUS**2.) \
                            * (planet_spectra.wavelength / h*c) * interp_function(planet_spectra.wavelength)
            else:
                signal      = (regridded_star_spectrum * STAR_RADIUS**2. + planet_spectra.flux * PLANET_RADIUS**2.) \
                            * (planet_spectra.wavelength / h*c) * 1.

                signal_star = (regridded_star_spectrum * STAR_RADIUS**2.) \
                            * (planet_spectra.wavelength / h*c) * 1.

            
            # Integrate the total signal and the star signal
            integrated_signal[i]      = trapz(signal, x=planet_spectra.wavelength)
            integrated_signal_star[i] = trapz(signal_star,  x=planet_spectra.wavelength)
            
        # Convert both the signal and the star signal to arrays 
        integrated_signal      = np.asarray(integrated_signal)
        integrated_signal_star = np.asarray(integrated_signal_star)
        
        # Divide out the two integrated signals
        fp_fs_ratio = (integrated_signal / integrated_signal_star - 1.0)
        
        # Save the data
        pd.DataFrame({'Phase':np.arange(0, 360, rot_val), 'fp_fs_ratio':fp_fs_ratio*1e6}).to_csv('OUTPUT_DATA/Phase_Curve_{}.txt'.format(planet_name), sep=' ')
        
        if ("SOOT".lower() in planet_name.lower()):
            linestyle_str = 'dotted'
        elif ("THOLIN".lower() in planet_name.lower()):
            linestyle_str = 'dashed'
        else:
            linestyle_str = 'solid'
        
        # Plot the data
        phases = np.arange(0, 360, 4*3.75) / 360
        ax.plot(phases,fp_fs_ratio * 1e6,
                linestyle=linestyle_str,
                color=my_colors(k / len(planet_names)),
                linewidth=3, label=planet_name)

    # Figure legend stuff
    ax.set_xlim(min(phases),max(phases))   
    ax.legend(fontsize=12, loc=(0,1.08), ncol=2, mode='expand')
    #ax.legend(fontsize=9, ncol=3, bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left", borderaxespad=0, mode="expand")
    ax.set_xlabel('Orbital Phase')
    ax.set_ylabel(r'F$_p$/F$_s$ (ppm)') #  (W m$^{-2}$)
    plt.savefig('../Figures/{}Phase_Curves.jpg'.format(planet_names[0][:planet_name_char_len]), dpi=200, bbox_inches='tight')
    plt.clf()
    return None
