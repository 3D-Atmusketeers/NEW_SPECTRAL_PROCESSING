#!/usr/bin/env python
# coding: utf-8

import pandas as pd

import re
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.pylab as pl
from scipy import interpolate
from scipy import ndimage
import numpy as np
import scipy.interpolate as spi
from scipy.integrate import trapz
from astropy import units as u
from astropy.modeling.models import BlackBody
from astropy.io import fits

# GLOBAL PARAMS
h = 6.6260701e-34
c = 3e8


#print ("Plotting the spectra...")
#print ()
#print ()

def filter_spectra(planet_spectra, star_spectra, filt, transmission_filter_name, wav_subset, phase_degrees):
    # Load in the planet spectra and convert to per microns
    #planet_spectra = pd.read_csv(file_path.format(str(phase_degrees)), header=None, delimiter=r"\s+", names=['wavelength', 'flux', 'reflected'])
    #planet_spectra.flux = planet_spectra.flux * (3.0e8 / planet_spectra.wavelength ** 2) / 1e6

    if transmission_filter_name != 'None':
        # Only use the data points within the filter wavelength range that we want
        planet_spectra = planet_spectra[(planet_spectra['wavelength'] > np.amin(filt.wav)) & (planet_spectra['wavelength'] < np.amax(filt.wav))]
        star_spectra = star_spectra[(star_spectra['wavelength'] > np.amin(filt.wav)) & (star_spectra['wavelength'] < np.amax(filt.wav))]

        # Cut off anything that is outside the wavelength range that you want
        planet_spectra = planet_spectra[(planet_spectra['wavelength'] > wav_subset[0]) & (planet_spectra['wavelength'] < wav_subset[1])]
        star_spectra = star_spectra[(star_spectra['wavelength'] > wav_subset[0]) & (star_spectra['wavelength'] < wav_subset[1])]

        # Reset the indexes of the planet and star dataframes
        planet_spectra = planet_spectra.reset_index(drop=True)
        star_spectra = star_spectra.reset_index(drop=True)

    return planet_spectra, star_spectra


def get_fp_fs(interp_function, star_radius, planet_spectra, planet_radius, star_spectra):
    # Re grid the star to the same wavelength grid as the planet
    fnc_star_converter = interpolate.interp1d(star_spectra.wavelength, star_spectra.flux, fill_value='extrapolate')
    regridded_star_spectrum = fnc_star_converter(planet_spectra.wavelength)

    # Ignore the filter if none is set
    if interp_function is None:
        signal = (regridded_star_spectrum * star_radius**2. + planet_spectra.flux * planet_radius ** 2.) * \
                 (planet_spectra.wavelength / h * c)
        signal_star = (regridded_star_spectrum * star_radius ** 2.) * (planet_spectra.wavelength / h * c)
    else:
        signal = (regridded_star_spectrum * star_radius**2. + planet_spectra.flux * planet_radius ** 2.) * \
                 (planet_spectra.wavelength / h * c) * interp_function(planet_spectra.wavelength)
        signal_star = (regridded_star_spectrum * star_radius ** 2.) * (planet_spectra.wavelength / h * c) * \
                       interp_function(planet_spectra.wavelength)


    return signal, signal_star


def get_filter(which_filter):
    if which_filter == 'MIRI':
        spectral_filter = pd.read_csv('DATA/MIRI_BANDPASS.txt', delimiter=r"\s+", skiprows=1, names=['wav', 'transmission'])
        spectral_filter.wav = spectral_filter.wav / 1e6
        spectral_filter_function = interpolate.interp1d(spectral_filter.wav, spectral_filter.transmission)
        return spectral_filter, spectral_filter_function
    elif which_filter == 'SPITZER_3_6':
        spectral_filter = pd.read_csv('DATA/filter_3_6.dat', delimiter=r"\s+", skiprows=1, names=['wav', 'transmission'])
        spectral_filter.wav = spectral_filter.wav / 1e6
        spectral_filter_function = interpolate.interp1d(spectral_filter.wav, spectral_filter.transmission)
        return spectral_filter, spectral_filter_function
    elif which_filter == 'SPITZER_4_5':
        spectral_filter = pd.read_csv('DATA/filter_4_5.dat', delimiter=r"\s+", skiprows=1, names=['wav', 'transmission'])
        spectral_filter.wav = spectral_filter.wav / 1e6
        spectral_filter_function = interpolate.interp1d(spectral_filter.wav, spectral_filter.transmission)
        return spectral_filter, spectral_filter_function
    else:
        print("No filter set!")
        return None, None


def get_star_spectra(planet_name):
    if "gj1214" in planet_name.replace("_", "").replace("-", "").lower():
        #print("Using a GJ1214 stellar spectrum")
        # Get the initial flux
        stellar_spectrum_1 = pd.read_csv(
            'DATA/GJ1214b_stellar_spectrum.txt',
            names=[
                'wavelength',
                'flux_Watt/m^2/micron'],
            delimiter=r"\s+",
            skiprows=1)
        wavelengths_meters = np.asarray(list(stellar_spectrum_1.wavelength * 1e-6)) * u.m
        flux_Watt_m_2_microns = np.asarray(list(stellar_spectrum_1['flux_Watt/m^2/micron'])) * u.W / u.m**2 / u.um
        flux_si = flux_Watt_m_2_microns * np.pi * u.sr
        star_spectra = pd.DataFrame(list(zip(wavelengths_meters.value, flux_si.value)), columns=['wavelength', 'flux'])

        #star_radius = 0.215 * 6.957e8
        star_radius = 1.50493470e8

    #elif ("peg51".lower() in planet_name.lower()):
    elif "peg51" in planet_name.replace("_", "").replace("-", "").lower():
        print("Using a peg 51b stellar spectra")
        spectra_file = 'DATA/peg51_lte05800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
        with fits.open(spectra_file) as hdul:
            spectra = hdul[0].data * 1e-7
        with fits.open('DATA/Phoenix_Wav.fits') as hdul:
            # Get the wavelengths in microns
            wavelengths = hdul[0].data * 1e-4

        stellar_spectrum_1 = pd.DataFrame(list(zip(wavelengths, spectra)), columns=['wavelength', 'flux_Watt/m^2/micron'])
        wavelengths_meters = np.asarray(list(stellar_spectrum_1.wavelength * 1e-6)) * u.m
        flux_Watt_m_2_microns = np.asarray(list(stellar_spectrum_1['flux_Watt/m^2/micron'])) * u.W / u.m**2 / u.um
        flux_si = flux_Watt_m_2_microns * u.sr
        star_spectra = pd.DataFrame(list(zip(wavelengths_meters.value, flux_si.value)), columns=['wavelength', 'flux'])

        star_radius = 1.237 * 6.957e8
    elif "taub" in planet_name.replace("_", "").replace("-", "").lower():
        print("Using a taub stellar spectra")
        spectra_file = 'DATA/taub_lte06400-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
        with fits.open(spectra_file) as hdul:
            spectra = hdul[0].data * 1e-7
        with fits.open('DATA/Phoenix_Wav.fits') as hdul:
            # Get the wavelengths in microns
            wavelengths = hdul[0].data * 1e-4

        stellar_spectrum_1 = pd.DataFrame(list(zip(wavelengths, spectra)), columns=['wavelength', 'flux_Watt/m^2/micron'])
        wavelengths_meters = np.asarray(list(stellar_spectrum_1.wavelength * 1e-6)) * u.m
        flux_Watt_m_2_microns = np.asarray(list(stellar_spectrum_1['flux_Watt/m^2/micron'])) * u.W / u.m**2 / u.um
        flux_si = flux_Watt_m_2_microns * u.sr
        star_spectra = pd.DataFrame(list(zip(wavelengths_meters.value, flux_si.value)), columns=['wavelength', 'flux'])
        star_radius = 1.42 * 6.957e8



    elif "wasp77a" in planet_name.replace("_", "").replace("-", "").lower():
        print("Using a wasp 77 a pheonix stellar spectra")
        spectra_file = 'DATA/lte05600-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
        with fits.open(spectra_file) as hdul:
            spectra = hdul[0].data * 1e-7
        with fits.open('DATA/Phoenix_Wav.fits') as hdul:
            # Get the wavelengths in microns
            wavelengths = hdul[0].data * 1e-4

        stellar_spectrum_1 = pd.DataFrame(list(zip(wavelengths, spectra)), columns=['wavelength', 'flux_Watt/m^2/micron'])
        wavelengths_meters = np.asarray(list(stellar_spectrum_1.wavelength * 1e-6)) * u.m
        flux_Watt_m_2_microns = np.asarray(list(stellar_spectrum_1['flux_Watt/m^2/micron'])) * u.W / u.m**2 / u.um
        flux_si = flux_Watt_m_2_microns * u.sr
        star_spectra = pd.DataFrame(list(zip(wavelengths_meters.value, flux_si.value)), columns=['wavelength', 'flux'])
        star_radius = 0.910 * 6.957e8


    elif "hd189" in planet_name.replace("_", "").replace("-", "").lower():
        print("Using a HD189 stellar spectrum")

        # Get the initial flux
        # This is in nm and ergs/s/cm**2/ster/nm
        stellar_spectrum = pd.read_csv('DATA/HD189733b_stellar_spectrum.txt',
                                       names=['wavelength', 'flux', 'hmm'],
                                       delimiter=r"\s+")
        wavelengths_meters = np.asarray(list(stellar_spectrum.wavelength * 1e-9)) * u.m

        # Interestingly, you get 1e-7 from ergs/s to Watts, and 1e4 * 1e3 from cm to m and nm to um
        flux_Watt_m_2_microns = np.asarray(list(stellar_spectrum['flux'])) * u.W / u.m**2 / u.um
        flux_si = flux_Watt_m_2_microns * u.sr * 4.0 * np.pi
        star_spectra = pd.DataFrame(list(zip(wavelengths_meters.value, flux_si.value)), columns=['wavelength', 'flux'])

        star_radius = 0.765 * 6.957e8

    elif "hd209" in planet_name.replace("_", "").replace("-", "").lower():
        print("Using a HD209 stellar spectrum")

        # Get the initial flux
        # This is in nm and ergs/s/cm**2/ster/nm
        stellar_spectrum = pd.read_csv(
            'DATA/HD209458b_stellar_spectrum.txt',
            names=['wavelength', 'flux', 'hmm'],
            delimiter=r"\s+")
        wavelengths_meters = np.asarray(list(stellar_spectrum.wavelength * 1e-9)) * u.m

        # Interestingly, you get 1e-7 from ergs/s to Watts, and 1e4 * 1e3 from cm to m and nm to um
        flux_Watt_m_2_microns = np.asarray(list(stellar_spectrum['flux'])) * u.W / u.m**2 / u.um
        flux_si = flux_Watt_m_2_microns * u.sr * 4.0 * np.pi
        star_spectra = pd.DataFrame(list(zip(wavelengths_meters.value, flux_si.value)), columns=['wavelength', 'flux'])

        star_radius = 1.203 * 6.957e8

    elif "wasp121" in planet_name.replace("_", "").replace("-", "").lower():
        print("Using a wasp121 stellar spectra")
        spectra_file = 'DATA/lte06500-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
        with fits.open(spectra_file) as hdul:
            spectra = hdul[0].data * 1e-7
        with fits.open('DATA/Phoenix_Wav.fits') as hdul:
            # Get the wavelengths in microns
            wavelengths = hdul[0].data * 1e-4

        stellar_spectrum_1 = pd.DataFrame(list(zip(wavelengths, spectra)), columns=['wavelength', 'flux_Watt/m^2/micron'])
        wavelengths_meters = np.asarray(list(stellar_spectrum_1.wavelength * 1e-6)) * u.m
        flux_Watt_m_2_microns = np.asarray(list(stellar_spectrum_1['flux_Watt/m^2/micron'])) * u.W / u.m**2 / u.um
        flux_si = flux_Watt_m_2_microns * u.sr
        star_spectra = pd.DataFrame(list(zip(wavelengths_meters.value, flux_si.value)), columns=['wavelength', 'flux'])

        star_radius = 0.805 * 6.957e8



    elif "hd189" in planet_name.replace("_", "").replace("-", "").lower():
        print("USING PHEONIX SPECTRA, BEWARE OF WAVELENGTH ENDS!!!!")
        spectra_file = 'DATA/lte04900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
        with fits.open(spectra_file) as hdul:
            spectra = hdul[0].data * 1e-7
        with fits.open('DATA/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits') as hdul:
            # Get the wavelengths in microns
            wavelengths = hdul[0].data * 1e-4

        stellar_spectrum_1 = pd.DataFrame(list(zip(wavelengths, spectra)), columns=['wavelength', 'flux_Watt/m^2/micron'])
        wavelengths_meters = np.asarray(list(stellar_spectrum_1.wavelength * 1e-6)) * u.m
        flux_Watt_m_2_microns = np.asarray(list(stellar_spectrum_1['flux_Watt/m^2/micron'])) * u.W / u.m**2 / u.um
        flux_si = flux_Watt_m_2_microns * u.sr
        star_spectra = pd.DataFrame(list(zip(wavelengths_meters.value, flux_si.value)), columns=['wavelength', 'flux'])
        star_radius = 1.458 * 6.957e8

    elif "55cnc" in planet_name.replace("_", "").replace("-", "").lower():
        print('USING A BLACKBODY FOR THE STAR')

        bb = BlackBody(temperature=5318 * u.K)
        wav = np.linspace(0.1, 20.0, 1000) * u.um
        flux = bb(wav) * u.sr * np.pi
        flux_si = flux.to(u.J / (u.m * u.m * u.Hz * u.s))
        wavelengths_meters = np.asarray(wav.value * 1e-6) * u.m
        star_spectra = pd.DataFrame(list(zip(wavelengths_meters.value, flux_si.value)), columns=['wavelength', 'flux'])
        star_spectra.flux = star_spectra.flux * (3.0e8 / star_spectra.wavelength ** 2) / 1e6
        star_radius = 0.96441 * 6.957e8


    else:
        print(planet_name)
        print("Your stellar spectrum isn't properly set")
        print("Hardwire here for a blackbody")
        exit(0)


    return star_spectra, star_radius


def specint(wn, spec, filtwn_list, filttrans_list):  # from Ryan Challener's Theresa code
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
    return (arr[-1] - arr[0]) / (len(arr) - 1)


def fwhm2sigma(fwhm):
    return fwhm / np.sqrt(8 * np.log(2))


def gaussian(x, mu, sig):
    normalization = 1 / np.sqrt(2 * np.pi * sig**2)
    exponent = - ((x - mu)**2 / (2 * sig**2))
    y = normalization * np.exp(exponent)
    return y


def reduceSpectralResolution(x, y, R_low, R_high=None, lambda_mid=None, n=4):
    if np.ndim(y) == 2:
        return np.array([reduceSpectralResolution(x, spec, R_low, R_high, lambda_mid, n) for spec in y])

    dx = getSpacing(x)
    # If lambda_mid is none, take median of input wavelengths
    if lambda_mid is None:
        lambda_mid = np.median(x)

    # If R_high is none, use midpoint of x divided by spacing in x
    if R_high is None:
        R_high = lambda_mid / dx

    # Create Gaussian kernel
    fwhm = np.sqrt(R_high**2 - R_low**2) * (lambda_mid / (R_low * R_high))
    sigma = fwhm2sigma(fwhm)

    kernel_x = np.arange(-n * sigma, n * sigma + dx, dx)
    kernel = gaussian(kernel_x, 0, sigma)
    kernel = kernel / np.sum(kernel)

    # find center of kernel
    n_kernel_lt0 = len(np.where(kernel_x < 0)[0])
    n_kernel_gt0 = len(np.where(kernel_x > 0)[0])

    if n_kernel_lt0 < n_kernel_gt0:
        origin = 0
    else:
        origin = -1

    # convolve
    lowRes = ndimage.convolve(y, kernel, origin=origin)
    return lowRes


# ## Spectra Test
def plot_planet_spectra_blackbody_comparison_hz(planet_names, black_body_temperatures, num_phases, INC_STRING):
    """
    plot the planet spectra, compared to a blackbody
    planet names: list of strings
    black_body_temperatures: list of floats
    """

    colors_bb = pl.cm.Greys(np.linspace(0, 1, len(black_body_temperatures) + 2))

    cm_name = 'batlow'
    cm_file = np.loadtxt(f'ScientificColourMaps8/{cm_name}/{cm_name}.txt')
    cm_file = np.roll(cm_file, 140, axis=0)
    my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file)


    for k, planet_name in enumerate(planet_names):
        # Figure aesthetics
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6), sharex=True, sharey=False)
        plt.subplots_adjust(hspace=0.05, wspace=0.25)


        for i in range(num_phases):
            # Get the phase value
            rot_val = 360. / num_phases

            # Get the file path
            file_path = '../Spectral-Processing/FINISHED_SPECTRA/Spec_0_' + planet_name + '_phase_{}_inc_' + INC_STRING + '00.00.dat'

            # Load in the planet spectra
            planet_spectra = pd.read_csv(file_path.format(str(i * rot_val)), header=None,
                                         delimiter=r"\s+", names=['wavelength', 'flux', 'reflected'])

            # Reset the index
            planet_spectra = planet_spectra.reset_index(drop=True)

            ax.plot(planet_spectra.wavelength * 1e6,
                    planet_spectra.flux,
                    alpha=1.0,
                    color=my_colors(i / num_phases),
                    linewidth=2.0,
                    label="Planet Spectra at phase: " + str(np.round(rot_val * i / 360., 3)))

        # Plot the blackbody temperatures
        j = 0
        for temp in black_body_temperatures:
            bb = BlackBody(temperature=temp * u.K)
            wav = np.linspace(0.5, 30.0, 1000) * u.um
            flux = bb(wav) * np.pi * u.sr

            # Convert to si, same as the output of the emission code
            flux = flux.to(u.J / (u.m * u.m * u.Hz * u.s))

            ax.plot(wav, flux,
                    color=colors_bb[j + 2],
                    alpha=1.0,
                    linewidth=2.0,
                    linestyle='dashed',
                    label='Blackbody at: ' + str(temp) + 'K')
            j = j + 1

        # Do some figure stuff
        #ax.set_ylim(1e-14, 5e-8)
        #ax.set_yscale('log')
        #ax.set_xlim(0.5, 5)
        ax.legend(fontsize=12, loc=(0, 1.05), ncol=2, mode='expand', title_fontsize=16)
        ax.set_xlabel(r'Wavelength ($\mu$m)')
        ax.set_ylabel(r'Flux (W/m$^2$Hz)')  # (W m$^{-2}$)
        plt.savefig('../Figures/planet_spectra_blackbody_comparison_hz_{}.jpg'.format(planet_name), dpi=200, bbox_inches='tight')
        plt.clf()
    return None



def plot_planet_spectra_blackbody_comparison_microns(planet_names, INC_STRING, black_body_temperatures, num_phases):
    """
    plot the planet spectra, compared to a blackbody
    planet names: list of strings
    black_body_temperatures: list of floats
    """

    colors_bb = pl.cm.Greys(np.linspace(0, 1, len(black_body_temperatures) + 2))

    cm_name = 'batlow'
    cm_file = np.loadtxt(f'ScientificColourMaps8/{cm_name}/{cm_name}.txt')
    cm_file = np.roll(cm_file, 140, axis=0)
    my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file)


    for k, planet_name in enumerate(planet_names):
        # Figure aesthetics
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6), sharex=True, sharey=False)
        plt.subplots_adjust(hspace=0.05, wspace=0.25)


        for i in range(num_phases):
            # Get the phase value
            rot_val = 360. / num_phases

            # Get the file path
            file_path = '../Spectral-Processing/FINISHED_SPECTRA/Spec_0_' + planet_name + '_phase_{}_inc_' + INC_STRING + '00.00.dat'

            # Load in the planet spectra
            planet_spectra = pd.read_csv(file_path.format(str(i * rot_val)), header=None,
                                         delimiter=r"\s+", names=['wavelength', 'flux', 'reflected'])

            # Reset the index
            planet_spectra = planet_spectra.reset_index(drop=True)

            # Convert to per micron
            planet_spectra.flux = planet_spectra.flux * (3.0e8 / planet_spectra.wavelength ** 2) / 1e6

            ax.plot(planet_spectra.wavelength * 1e6,
                    planet_spectra.flux,
                    alpha=1.0,
                    color=my_colors(i / num_phases),
                    linewidth=0.5,
                    label="Planet Spectra at phase: " + str(np.round(rot_val * i / 360., 3)))

        # Plot the blackbody temperatures
        j = 0
        for temp in black_body_temperatures:
            # Set up the blackbody at a given temperature
            bb = BlackBody(temperature=temp * u.K)

            # The wavelengths to use
            wav = np.linspace(0.5, 30.0, 1000) * u.um

            # Get the flux at that wavelength
            flux = bb(wav) * np.pi * u.sr

            # Convert to si, same as the output of the emission code
            flux = flux.to(u.J / (u.m * u.m * u.Hz * u.s))

            # Convert the wavelengths to meters
            wav = wav * 1e-6

            # Create a dataframe called bb_df out of the lists wav and flux
            bb_df = pd.DataFrame({'wavelength': wav, 'flux': flux})
            bb_df.flux = bb_df.flux * (3.0e8 / bb_df.wavelength ** 2) / 1e6

            ax.plot(bb_df.wavelength * 1e6, bb_df.flux,
                    color=colors_bb[j + 2],
                    alpha=1.0,
                    linewidth=2.0,
                    linestyle='dashed',
                    label='Blackbody at: ' + str(temp) + 'K')
            j = j + 1

        # Do some figure stuff
        #ax.set_ylim(1e3, 1e5)
        #ax.set_ylim(0, 1200)
        #ax.set_xlim(2, 15)
        #ax.set_yscale('log')
        ax.legend(fontsize=12, loc=(0, 1.05), ncol=2, mode='expand', title_fontsize=16)
        ax.set_xlabel(r'Wavelength ($\mu$m)')
        ax.set_ylabel(r'Flux (W/m$^2$micron)')  # (W m$^{-2}$)
        plt.savefig('../Figures/planet_spectra_blackbody_comparison_microns{}.jpg'.format(planet_name), dpi=200, bbox_inches='tight')
        plt.clf()
    return None


def plot_star_spectra_test(planet_names):
    # Figure aesthetics
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6), sharex=True, sharey=False)
    plt.subplots_adjust(hspace=0.05, wspace=0.25)

    for planet_name in planet_names:
        # get the star_spectra
        star_spectra, star_radius = get_star_spectra(planet_name)

        plt.plot(
            star_spectra.wavelength * 1e6,
            star_spectra.flux,
            label='Stellar Spectrum',
            alpha=1.0,
            color='black',
            linewidth=0.5)

        # blackbody comparison
        Teffs = [6500]
        for Teff in Teffs:
            bb = BlackBody(temperature=Teff * u.K)
            wav = np.linspace(0.1, 20.0, 1000) * u.um
            flux = bb(wav) * u.sr * np.pi
            flux_si = flux.to(u.J / (u.m * u.m * u.Hz * u.s))
            wavelengths_meters = np.asarray(wav.value * 1e-6) * u.m
            bb_star_spectra = pd.DataFrame(list(zip(wavelengths_meters.value, flux_si.value)), columns=['wavelength', 'flux'])

            # Convert it to microns
            bb_star_spectra.flux = bb_star_spectra.flux * (3.0e8 / bb_star_spectra.wavelength ** 2) / 1e6


            plt.plot(
                bb_star_spectra.wavelength *
                1e6,
                bb_star_spectra.flux,
                label='Blackbody at ' +
                str(Teff),
                alpha=1.0,
                color='red',
                linewidth=1.5)

        plt.ylim(bottom=1e2)
        plt.xlim(left=0.0)
        plt.yscale('log')
        plt.xlabel('Wavelengths (microns)')
        plt.ylabel(r'Flux (Watts/m$^2$/micron)')
        plt.legend(loc='upper right')
        plt.savefig('../Figures/star_test_{}.jpg'.format(planet_name), dpi=200, bbox_inches='tight')
        plt.clf()

    return None


def plot_filters(planet_names, transmission_filter_name, INC_STRING):
    spectral_filter, spectral_filter_function = get_filter(which_filter=transmission_filter_name)

    # Figure aesthetics
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
    plt.subplots_adjust(hspace=0.05, wspace=0.25)

    # Get the file path

    file_path = '../Spectral-Processing/FINISHED_SPECTRA/Spec_0_' + planet_names[0] + '_phase_{}_inc_' + INC_STRING + '00.00.dat'
    spectra = pd.read_csv(
        file_path.format(
            str(0.0)),
        header=None,
        delimiter=r"\s+",
        names=[
            'wavelength',
            'flux',
            'reflected'])

    filtered_spectra = spectra[(spectra['wavelength'] > np.amin(spectral_filter.wav)) & (spectra['wavelength'] < np.amax(spectral_filter.wav))]

    ax.plot(spectral_filter.wav * 1e6, spectral_filter.transmission, color='red', linewidth=3, alpha=0.9, label=transmission_filter_name)
    ax.plot(filtered_spectra.wavelength * 1e6, spectral_filter_function(filtered_spectra.wavelength),
            color='black',
            linewidth=3,
            linestyle='dashed',
            label="Interpolated")
        
    # Some more Figure stuff
    ax.set_xlabel(r"Wavelength ($\mu$m)")
    ax.set_ylabel('Transmittance')
    ax.legend(fontsize=14)
    plt.savefig('../Figures/Filters_{}.jpg'.format(transmission_filter_name), dpi=100, bbox_inches='tight')
    #return None


def plot_spectra_simple(planet_names, num_phases, INC_STRING):
    cm_name = 'batlow'
    cm_file = np.loadtxt(f'ScientificColourMaps8/{cm_name}/{cm_name}.txt')
    cm_file = np.roll(cm_file, 140, axis=0)
    my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file)

    for k, planet_name in enumerate(planet_names):
        # Figure aesthetics
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6), sharex=True, sharey=False)
        plt.subplots_adjust(hspace=0.05, wspace=0.25)

        for i in range(num_phases):
            # Get the phase value
            rot_val = 360. / num_phases

            # Get the file path

            file_path = '../Spectral-Processing/FINISHED_SPECTRA/Spec_0_' + planet_name + '_phase_{}_inc_' + INC_STRING + '00.00.dat'

            # Load in the planet spectra
            planet_spectra = pd.read_csv(file_path.format(str(i * rot_val)), header=None,
                                         delimiter=r"\s+", names=['wavelength', 'flux', 'reflected'])
            planet_spectra.flux = planet_spectra.flux * (3.0e8 / planet_spectra.wavelength ** 2) / 1e6

            # Reset the index
            planet_spectra = planet_spectra.reset_index(drop=True)




            # Replace the substring DOGRAY in planetnames with PICKET
            planet_name2 = planet_name.replace("DOGRAY", "PICKET")
            file_path2 = '../Spectral-Processing/FINISHED_SPECTRA/Spec_0_' + planet_name2 + '_phase_{}_inc_' + INC_STRING + '00.00.dat'

            # Load in the planet spectra
            planet_spectra2 = pd.read_csv(file_path2.format(str(i * rot_val)), header=None,
                                         delimiter=r"\s+", names=['wavelength', 'flux', 'reflected'])
            planet_spectra2.flux = planet_spectra2.flux * (3.0e8 / planet_spectra2.wavelength ** 2) / 1e6

            # Reset the index
            planet_spectra2 = planet_spectra2.reset_index(drop=True)

            ax.plot(planet_spectra.wavelength * 1e6,
                    planet_spectra.flux,
                    color=my_colors(i / num_phases),
                    linewidth=1,
                    label=str(np.round(rot_val * i / 360., 3)))

        # Do somet figure stuff
        #ax.set_xlim(min(planet_spectra.wavelength * 1e6), max(planet_spectra.wavelength * 1e6))
        ax.set_xlim(0.5, 5)
        ax.set_ylim(1e1,5e6)
        ax.legend(fontsize=12, loc=(0, 1.03), ncol=5, mode='expand', title='Orbital Phase', title_fontsize=16)
        ax.set_xlabel(r'Wavelength ($\mu$m)')
        ax.set_ylabel(r'Flux (W/m$^2$/micron)')
        #ax.set_yscale('log')
        plt.savefig('../Figures/Planet_Simple_Spectra_{}.jpg'.format(planet_name), dpi=200, bbox_inches='tight')
        plt.clf()

    return None


def plot_fp_fs_phase_curves(planet_names, planet_name_char_len, planet_radii, num_phases,
                            transmission_filter_name, wav_subset, INC_STRING):
    # Figure aesthetics
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6), sharex=True, sharey=False)
    plt.subplots_adjust(hspace=0.05, wspace=0.25)

    cm_name = 'batlow'
    cm_file = np.loadtxt(f'ScientificColourMaps8/{cm_name}/{cm_name}.txt')
    cm_file = np.roll(cm_file, 0, axis=0)
    my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file)

    # I currently only have MIRI coded in
    filt, interp_function = get_filter(which_filter=transmission_filter_name)

    for k, planet_name in enumerate(planet_names):
        planet_radius = planet_radii[k]

        # Set the arrays for the integrated signals
        integrated_signal = []
        integrated_signal_star = []

        # Get the star spectra
        star_spectra, star_radius = get_star_spectra(planet_name)

        print('Setting the wavelength subset for the filter')

        # Extract wavelengths for ease
        start_star = list(star_spectra.wavelength)[0]
        end_star = list(star_spectra.wavelength)[-1]

        # Check against star_spectra and update wav_subset
        if start_star > wav_subset[0]:
            wav_subset[0] = start_star
            print("Adjusted the lower bound of wavelength subset based on star_spectra")

        if end_star < wav_subset[1]:
            wav_subset[1] = end_star
            print("Adjusted the upper bound of wavelength subset based on star_spectra")

        # If condition is met, further adjust wav_subset by comparing against filt data
        if transmission_filter_name.lower() not in ['none']:
            start_filt, end_filt = list(filt.wav)[0], list(filt.wav)[-1]
            wav_subset[0] = max(start_filt, wav_subset[0])
            wav_subset[1] = min(end_filt, wav_subset[1])

        for i in range(num_phases):
            # For reading in the file names
            rot_val = 360. / num_phases
            phase_degrees = rot_val * i

            # Add a zero as a placeholder
            integrated_signal.append(0)
            integrated_signal_star.append(0)

            # Get the file path
            file_path = '../Spectral-Processing/FINISHED_SPECTRA/Spec_0_' + planet_name + '_phase_{}_inc_' + INC_STRING + '00.00.dat'

            # Load in the planet spectra
            planet_spectra = pd.read_csv(file_path.format(str(i * rot_val)), header=None,
                                         delimiter=r"\s+", names=['wavelength', 'flux', 'reflected'])
            planet_spectra.flux = planet_spectra.flux * (3.0e8 / planet_spectra.wavelength ** 2) / 1e6

            # Reset the index
            planet_spectra = planet_spectra.reset_index(drop=True)

            # Call the function to get the filtered star and planet spectra
            planet_spectra, star_spectra = filter_spectra(planet_spectra, star_spectra, filt, transmission_filter_name, wav_subset, phase_degrees)

            # Get the star and planet plus star signals
            signal, signal_star = get_fp_fs(interp_function, star_radius, planet_spectra, planet_radius, star_spectra)

            # Integrate the total signal and the star signal
            integrated_signal[i] = trapz(signal, x=planet_spectra.wavelength * 1e6)
            integrated_signal_star[i] = trapz(signal_star, x=planet_spectra.wavelength * 1e6)

        # Convert both the signal and the star signal to arrays
        integrated_signal = np.asarray(integrated_signal)
        integrated_signal_star = np.asarray(integrated_signal_star)

        # Divide out the two integrated signals
        fp_fs_ratio = (integrated_signal / integrated_signal_star - 1.0)

        # Save the data
        pd.DataFrame({'Phase': np.arange(0, 360, rot_val), 'Fp_Fs_pmm': fp_fs_ratio * 1e6}
                     ).to_csv('OUTPUT_DATA/Fp_Fs_{}_Phase_Curves_{}.txt'.format(planet_name, transmission_filter_name), sep=' ')

        if '50clouds' in planet_name.lower():
            linestyle_str = 'dashed'
        elif '0clouds' in planet_name.lower():
            linestyle_str = 'solid'
        else:
            linestyle_str = 'dashdot'

        if k % 3 == 0:
            color_str='#70A5FF'
        elif k % 3 == 1:
            color_str='#832071'
        else:
            color_str='#E58B9D'

        # Plot the data
        phases = np.linspace(0, 345, num_phases) / 360
        ax.plot(phases, fp_fs_ratio * 1e6,
                #color='black',
                #linestyle=linestyle_str,
                linewidth=2,
                label=planet_name)

    
    """
    # Load the data
    df = pd.read_csv('/home/imalsky/Documents/Spectra-Paper/DATA/HD189bo11.csv')

    # Subtract 1 from 'phase' values that are greater than 1
    df['phase'] = df['phase'].apply(lambda x: x - 1 if x > 1 else x)

    # Prepare y_data
    y_data = (df.bestfit_norm[::100] - 1) * 1e6
    
    # Create the scatter plot
    ax.scatter(df.phase[::100], y_data, s=1)
    ax.set_ylim(0, 2500)
    """

    #ax.set_xlim(0, 0.15)
    #ax.set_ylim(0.1, 2000)

    ax.set_xlim(min(phases), max(phases))
    ax.legend(fontsize=9, loc=(0, 1.03), ncol=3, mode='expand')
    ax.set_xlabel('Orbital Phase')
    ax.set_ylabel(r'F$_p$/F$_s$ (ppm)')
    plt.savefig('../Figures/Fp_Fs_Phase_Curves_{}.jpg'.format(transmission_filter_name), dpi=200, bbox_inches='tight')
    plt.clf()
    return None


def plot_fp_phase_curves(planet_names, planet_name_char_len, num_phases,
                         transmission_filter_name, wav_subset, INC_STRING):
    # Figure aesthetics
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6), sharex=True, sharey=False)
    plt.subplots_adjust(hspace=0.05, wspace=0.25)

    cm_name = 'batlow'
    cm_file = np.loadtxt(f'ScientificColourMaps8/{cm_name}/{cm_name}.txt')
    cm_file = np.roll(cm_file, 140, axis=0)
    my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file)

    # I currently only have MIRI coded in
    filt, interp_function = get_filter(which_filter=transmission_filter_name)

    # Set the arrays for the integrated signals
    # Based on the filter
    print('Setting the wavelength subset for the filter')
    if (transmission_filter_name != 'None' and transmission_filter_name != 'none'):
        wav_subset[0] = list(filt.wav)[0]
        wav_subset[1] = list(filt.wav)[-1]
    else:
        pass

    for k, planet_name in enumerate(planet_names):
        # Set the arrays for the integrated signals
        integrated_planet_spectra = []

        # Get the star spectra
        star_spectra, star_radius = get_star_spectra(planet_name)

        running_flux_sum = []
        for i in range(num_phases):
            # For reading in the file names
            rot_val = 360. / num_phases
            phase_degrees = rot_val * i

            # Add a zero as a placeholder
            integrated_planet_spectra.append(0)

            # Load in the planet spectra and convert to per microns
            # Get the file path
            file_path = '../Spectral-Processing/FINISHED_SPECTRA/Spec_0_' + planet_name + '_phase_{}_inc_' + INC_STRING + '00.00.dat'

            # Load in the planet spectra
            planet_spectra = pd.read_csv(file_path.format(str(i * rot_val)), header=None,
                                         delimiter=r"\s+", names=['wavelength', 'flux', 'reflected'])
            
            planet_spectra['flux_all']   = planet_spectra.flux      * (3.0e8 / planet_spectra.wavelength ** 2) / 1e6
            planet_spectra['flux_refl']  = planet_spectra.reflected * (3.0e8 / planet_spectra.wavelength ** 2) / 1e6
            planet_spectra['flux']       = planet_spectra.flux_all - planet_spectra.flux_refl

            running_flux_sum.append(np.sum(planet_spectra.flux))
            

            # Reset the index
            planet_spectra = planet_spectra.reset_index(drop=True)

            # Call the function to get the filtered star and planet spectra
            planet_spectra, star_spectra = filter_spectra(planet_spectra, star_spectra, filt, transmission_filter_name, wav_subset, phase_degrees)


            if interp_function is None:
                planet_spectra.flux = planet_spectra.flux
            else:
                planet_spectra.flux = planet_spectra.flux * interp_function(planet_spectra.wavelength)

            # Integrate the total planet planet flux
            integrated_planet_spectra[i] = trapz(planet_spectra.flux, x=planet_spectra.wavelength * 1e6)
        
        #running_flux_sum = np.asarray(running_flux_sum)
        #print('Flux running sum', np.sum(running_flux_sum))

        # Convert both the signal and the star signal to arrays
        integrated_planet_spectra = np.asarray(integrated_planet_spectra)

        # Divide out the two integrated signals
        fp = integrated_planet_spectra

        # Save the data
        pd.DataFrame({'Phase': np.arange(0, 360, rot_val), 'fp W/m2': fp}
                     ).to_csv('OUTPUT_DATA/Fp_{}_Phase_Curves.txt'.format(planet_name), sep=' ')

        # Plot the data
        phases = np.linspace(0, 345, num_phases) / 360
        ax.plot(phases, fp,
                linestyle='solid',
                color=my_colors(k / len(planet_names)),
                linewidth=2.0,
                label=planet_name)

    # Figure legend stuff
    ax.set_xlim(0, max(phases))
    ax.set_ylim(0, 1e4)
    ax.legend(fontsize=12, loc=(0, 1.03), ncol=2, mode='expand')
    ax.set_xlabel('Orbital Phase')
    ax.set_ylabel(r'Planet Flux (W/m$^2$)')
    plt.savefig('../Figures/Fp_Phase_Curves_{}.jpg'.format(transmission_filter_name), dpi=200, bbox_inches='tight')
    plt.clf()
    return None





def plot_dayside(planet_names, planet_radii, num_phases, transmission_filter_name,
                       wav_subset, resolution, INC_STRING):
    cm_name = 'batlow'
    cm_file = np.loadtxt(f'ScientificColourMaps8/{cm_name}/{cm_name}.txt')
    cm_file = np.roll(cm_file, 140, axis=0)
    my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file)

    # I currently only have MIRI coded in
    filt, interp_function = get_filter(which_filter=transmission_filter_name)

    # Figure aesthetics
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6), sharex=True, sharey=False)
    plt.subplots_adjust(hspace=0.05, wspace=0.25)

    for k, planet_name in enumerate(planet_names):
        planet_radius = planet_radii[k]

        # Get the star spectra
        star_spectra, star_radius = get_star_spectra(planet_name)

        # Extract wavelengths for ease
        start_star = list(star_spectra.wavelength)[0]
        end_star = list(star_spectra.wavelength)[-1]

        # Check against star_spectra and update wav_subset
        if start_star > wav_subset[0]:
            wav_subset[0] = start_star
            print("Adjusted the lower bound of wavelength subset based on star_spectra")

        if end_star < wav_subset[1]:
            wav_subset[1] = end_star
            print("Adjusted the upper bound of wavelength subset based on star_spectra")

        # If condition is met, further adjust wav_subset by comparing against filt data
        if transmission_filter_name.lower() not in ['none']:
            start_filt, end_filt = list(filt.wav)[0], list(filt.wav)[-1]
            wav_subset[0] = max(start_filt, wav_subset[0])
            wav_subset[1] = min(end_filt, wav_subset[1])

        phase_degrees = 180.0

        # Get the file path
        file_path = '../Spectral-Processing/FINISHED_SPECTRA/Spec_0_' + planet_name + '_phase_{}_inc_' + INC_STRING + '00.00.dat'

        # Load in the planet spectra
        planet_spectra = pd.read_csv(file_path.format(str(180.0)), header=None,
                                        delimiter=r"\s+", names=['wavelength', 'flux', 'reflected'])
        planet_spectra.flux = planet_spectra.flux * (3.0e8 / planet_spectra.wavelength ** 2) / 1e6

        # Reset the index
        planet_spectra = planet_spectra.reset_index(drop=True)

        # Call the function to get the filtered star and planet spectra
        planet_spectra, star_spectra = filter_spectra(planet_spectra, star_spectra, filt, transmission_filter_name, wav_subset, phase_degrees)

        # Get the planet to star flux ratio
        signal, signal_star = get_fp_fs(interp_function, star_radius, planet_spectra, planet_radius, star_spectra)

        # Divide out the two integrated signals
        fp_fs_ratio = (signal / signal_star - 1.0)

        if resolution > 0:
            fp_fs_ratio = np.asarray(reduceSpectralResolution(list(planet_spectra.wavelength),list(fp_fs_ratio),
                                                                R_low=resolution, R_high=None, lambda_mid=None, n=4))
            
        pd.DataFrame({'Wavelength (microns)': planet_spectra.wavelength * 1e6, 'Fp_Fs_pmm': fp_fs_ratio * 1e2}
                         ).to_csv('OUTPUT_DATA/Fp_Fs_Spectra_{}_Spectra_{}.txt'.format(str(180.0), planet_name), sep=' ')

        ax.plot(planet_spectra.wavelength * 1e6,
                fp_fs_ratio * 1e2,
                #color=my_colors(k / len(planet_names)),
                linewidth=2,
                zorder=1,
                label=planet_name)
        


    df = pd.read_csv('/home/marianne/Desktop/Spectra-Paper/DATA/table_HD-189733-b-Changeat-et-al.-2022.csv',
                        skiprows=1, names=['wav', 'band', 'depth', 'err1', 'err2',
                                            'blank1', 'blank2', 'blank3', 'blank4'])

    # Set up y-error as a tuple of absolute values for lower and upper errors
    yerr = (df['err2'].abs(), df['err1'].abs())

    # Create the error bar plot with both y- and x-error
    ax.errorbar(
        df['wav'],
        df['depth'],
        #xerr=xerr,
        yerr=yerr,
        fmt='o',
        linewidth=1,
        color='blue',
        label='Changeat (2022)')

    # Figure legend
    ax.minorticks_on()  # This enables minor ticks
    ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.set_ylim(bottom=3e-3)
    ax.legend(fontsize=12, loc=(0, 1.03), ncol=2, mode='expand', title='Secondary Eclipse Emission Spectra', title_fontsize=12)
    ax.set_xlabel(r'Wavelength ($\mu$m)')
    ax.set_ylabel(r'F$_p$/F$_s$ (%)')  # (W m$^{-2}$)
    plt.savefig('../Figures/TEST.jpg', dpi=200, bbox_inches='tight')













def plot_fp_fs_spectra(planet_names, planet_radii, num_phases, transmission_filter_name,
                       wav_subset, resolution, INC_STRING):
    cm_name = 'batlow'
    cm_file = np.loadtxt(f'ScientificColourMaps8/{cm_name}/{cm_name}.txt')
    cm_file = np.roll(cm_file, 140, axis=0)
    my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file)

    # I currently only have MIRI coded in
    filt, interp_function = get_filter(which_filter=transmission_filter_name)

    for k, planet_name in enumerate(planet_names):
        # Figure aesthetics
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6), sharex=True, sharey=False)
        plt.subplots_adjust(hspace=0.05, wspace=0.25)

        planet_radius = planet_radii[k]

        # Get the star spectra
        star_spectra, star_radius = get_star_spectra(planet_name)

        print('Setting the wavelength subset for the filter')

        # Extract wavelengths for ease
        start_star = list(star_spectra.wavelength)[0]
        end_star = list(star_spectra.wavelength)[-1]

        # Check against star_spectra and update wav_subset
        if start_star > wav_subset[0]:
            wav_subset[0] = start_star
            print("Adjusted the lower bound of wavelength subset based on star_spectra")

        if end_star < wav_subset[1]:
            wav_subset[1] = end_star
            print("Adjusted the upper bound of wavelength subset based on star_spectra")

        # If condition is met, further adjust wav_subset by comparing against filt data
        if transmission_filter_name.lower() not in ['none']:
            start_filt, end_filt = list(filt.wav)[0], list(filt.wav)[-1]
            wav_subset[0] = max(start_filt, wav_subset[0])
            wav_subset[1] = min(end_filt, wav_subset[1])

        for i in range(num_phases):
            rot_val = 360. / num_phases
            phase_degrees = rot_val * i

            # Get the file path
            file_path = '../Spectral-Processing/FINISHED_SPECTRA/Spec_0_' + planet_name + '_phase_{}_inc_' + INC_STRING + '00.00.dat'

            # Load in the planet spectra
            planet_spectra = pd.read_csv(file_path.format(str(i * rot_val)), header=None,
                                         delimiter=r"\s+", names=['wavelength', 'flux', 'reflected'])
            planet_spectra.flux = planet_spectra.flux * (3.0e8 / planet_spectra.wavelength ** 2) / 1e6


            # Reset the index
            planet_spectra = planet_spectra.reset_index(drop=True)

            # Call the function to get the filtered star and planet spectra
            planet_spectra, star_spectra = filter_spectra(planet_spectra, star_spectra, filt, transmission_filter_name, wav_subset, phase_degrees)
            
            # Get the planet to star flux ratio
            signal, signal_star = get_fp_fs(interp_function, star_radius, planet_spectra, planet_radius, star_spectra)
            #signal, signal_star = get_fp_fs(interp_function, star_radius, blackbody, planet_radius, star_spectra)

            # Divide out the two integrated signals
            fp_fs_ratio = (signal / signal_star - 1.0)

            if resolution > 0:
                fp_fs_ratio = np.asarray(reduceSpectralResolution(list(planet_spectra.wavelength),list(fp_fs_ratio),
                                                                  R_low=resolution, R_high=None, lambda_mid=None, n=4))


            ax.plot(planet_spectra.wavelength.values * 1e6,
                    fp_fs_ratio * 1e6,
                    color=my_colors(i / num_phases),
                    linewidth=2,
                    zorder=1,
                    label=str(np.round(rot_val * i / 360., 3)))

            pd.DataFrame({'Wavelength (microns)': planet_spectra.wavelength * 1e6, 'Fp_Fs_pmm': fp_fs_ratio * 1e6}
                         ).to_csv('OUTPUT_DATA/Fp_Fs_Spectra_{}_Spectra_{}.txt'.format(str(i * rot_val), planet_name), sep=' ')
            
            #pd.DataFrame({'Wavelength (microns)': planet_spectra.wavelength * 1e6, 'Fp_Fs_pmm': fp_fs_ratio * 1e6}
            #             ).to_csv('OUTPUT_DATA/Blackbody_test.txt', sep=' ')


            #spectrum_data = pd.read_csv('/home/marianne/Desktop/HD189733b_data.txt',
            #                            delimiter=r"\s+",
            #                            skiprows=1,
            #                            names=['Wavelength (microns)', 'Fp_Fs_pmm', 'error'])
            #ax.errorbar(
            #    spectrum_data['Wavelength (microns)'],  # X-axis
            #    spectrum_data['Fp_Fs_pmm'],             # Y-axis
            ##    yerr=spectrum_data['error'],            # Error values
            #    fmt='o',                               # Format (circle markers with line)
            #    linewidth=1,                            # Line width
            #    label='Observational Data',  # Label for the plot
            #    zorder=2)

        
        # Figure legend
        #ax.set_ylim(1, 1399)
        ax.minorticks_on()  # This enables minor ticks
        ax.set_xlim(min(planet_spectra.wavelength * 1e6), max(planet_spectra.wavelength * 1e6))
        #ax.set_yscale('log')
        #ax.set_xscale('log')
        #ax.set_xlim(5,13)
        #ax.set_ylim(1e2, 1e6)
        #ax.set_ylim(bottom=1e-2)
        ax.legend(fontsize=12, loc=(0, 1.03), ncol=5, mode='expand', title='Orbital Phase', title_fontsize=18)
        ax.set_xlabel(r'Wavelength ($\mu$m)')
        ax.set_ylabel(r'F$_p$/F$_s$ (ppm)')  # (W m$^{-2}$)
        plt.savefig('../Figures/Fp_Fs_Spectra_{}_{}.jpg'.format(planet_name, transmission_filter_name), dpi=200, bbox_inches='tight')

    return None


def plot_fp_spectra(planet_names, num_phases, transmission_filter_name, wav_subset, resolution, INC_STRING):
    cm_name = 'batlow'
    cm_file = np.loadtxt(f'ScientificColourMaps8/{cm_name}/{cm_name}.txt')
    cm_file = np.roll(cm_file, 140, axis=0)
    my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file)

    # I currently only have MIRI coded in
    filt, interp_function = get_filter(which_filter=transmission_filter_name)

    # Set the arrays for the integrated signals
    # Based on the filter
    print('Setting the wavelength subset for the filter')
    if (transmission_filter_name != 'None' and transmission_filter_name != 'none'):
        wav_subset[0] = list(filt.wav)[0]
        wav_subset[1] = list(filt.wav)[-1]
    else:
        pass

    for k, planet_name in enumerate(planet_names):
        # Figure aesthetics
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6), sharex=True, sharey=False)
        plt.subplots_adjust(hspace=0.05, wspace=0.25)

        # Get the star spectra
        star_spectra, star_radius = get_star_spectra(planet_name)

        for i in range(num_phases):
            rot_val = 360. / num_phases
            phase_degrees = rot_val * i

            # Get the file path
            file_path = '../Spectral-Processing/FINISHED_SPECTRA/Spec_0_' + planet_name + '_phase_{}_inc_' + INC_STRING + '00.00.dat'

            # Load in the planet spectra
            planet_spectra = pd.read_csv(file_path.format(str(i * rot_val)), header=None,
                                         delimiter=r"\s+", names=['wavelength', 'flux', 'reflected'])
                        
            planet_spectra.flux = planet_spectra.flux * (3.0e8 / planet_spectra.wavelength ** 2) / 1e6


            # Reset the index
            planet_spectra = planet_spectra.reset_index(drop=True)

            # Call the function to get the filtered star and planet spectra
            planet_spectra, star_spectra = filter_spectra(planet_spectra, star_spectra, filt, transmission_filter_name, wav_subset, phase_degrees)

            if interp_function is None:
                planet_spectra.flux = planet_spectra.flux
            else:
                planet_spectra.flux = planet_spectra.flux * interp_function(planet_spectra.wavelength)

            flux = planet_spectra.flux

            if resolution > 0:
                flux = np.asarray(reduceSpectralResolution(list(planet_spectra.wavelength),list(flux),
                                                           R_low=resolution, R_high=None, lambda_mid=None, n=4))

            ax.plot(planet_spectra.wavelength.values * 1e6,
                    list(flux),
                    color=my_colors(i / num_phases),
                    linewidth=2.0,
                    label=str(np.round(rot_val * i / 360., 3)))

            pd.DataFrame({'Wavelength (microns)': planet_spectra.wavelength * 1e6, 'Planet_Flux_micron': flux}
                         ).to_csv('OUTPUT_DATA/Fp_Spectra_{}_{}_{}.txt'.format(str(i * rot_val), planet_name, transmission_filter_name), sep=' ')


        """
        df = pd.read_excel('/home/imalsky/Desktop/41586_2023_6159_MOESM2_ESM.xlsx',
                           skiprows=3,
                           names=["Wavelength","Fp_day","Up 1-sigma","Low 1-sigma","Fp_night","Up 1-sigma","Low 1-sigma"])

        ax.errorbar(df['Wavelength'], df['Fp_day'], yerr=[df['Low 1-sigma'], df['Up 1-sigma']],
             fmt='o-', label='Fp_day', capsize=2, linewidth=2)
        ax.errorbar(df['Wavelength'], df['Fp_night'], yerr=[df['Low 1-sigma.1'], df['Up 1-sigma.1']],
             fmt='o-', label='Fp_night', capsize=2, linewidth=2)
        """
        #ax.set_xscale('log')
        #ax.set_yscale('log') 
        #ax.set_xlim(2.3105, 2.3135)
        #ax.set_ylim(2e2, 2e5)
        ax.minorticks_on()  # This enables minor ticks

        legend = ax.legend(fontsize=12, loc=(0, 1.03), ncol=6, title='Orbital Phase', mode='expand',title_fontsize=18)
        for line in legend.get_lines():
            line.set_linewidth(3)  # Set the desired linewidth here

        ax.set_xlabel(r'Wavelength ($\mu$m)')
        ax.set_ylabel(r'F$_p$ (W/m$^2$/micron)')  # (W m$^{-2}$)
        plt.savefig('../Figures/Fp_Spectra_{}_{}.jpg'.format(planet_name, transmission_filter_name), dpi=200, bbox_inches='tight')


    return None


def plot_blackbody_phase_curve(planet_name, planet_radii,num_phases,transmission_filter_name, wav_subset,resolution,temp):
    cm_name = 'batlow'
    cm_file = np.loadtxt(f'ScientificColourMaps8/{cm_name}/{cm_name}.txt')
    cm_file = np.roll(cm_file, 140, axis=0)
    my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file)

    # I currently only have MIRI coded in
    filt, interp_function = get_filter(which_filter=transmission_filter_name)

    # Figure aesthetics
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6), sharex=True, sharey=False)
    plt.subplots_adjust(hspace=0.05, wspace=0.25)

    # Set up the blackbody at a given temperature
    bb = BlackBody(temperature=temp * u.K)

    # The wavelengths to use
    wav = np.linspace(5, 12, 1000) * u.um

    # Get the flux at that wavelength
    flux = bb(wav) * np.pi * u.sr

    # Convert to si, same as the output of the emission code
    flux = flux.to(u.J / (u.m * u.m * u.Hz * u.s))

    # Convert the wavelengths to meters
    wav = wav * 1e-6

    # Create a dataframe called bb_df out of the lists wav and flux
    bb_spectra = pd.DataFrame({'wavelength': wav, 'flux': flux})
    bb_spectra.flux = bb_spectra.flux * (3.0e8 / bb_spectra.wavelength ** 2) / 1e6

    planet_radius = planet_radii[0]

    # Get the star spectra
    star_spectra, star_radius = get_star_spectra(planet_name)

    # Set the arrays for the integrated signals
    integrated_signal = []
    integrated_signal_star = []

    # Get the star spectra
    star_spectra, star_radius = get_star_spectra(planet_name)

    for i in range(num_phases):
        # For reading in the file names
        rot_val = 360. / num_phases
        phase_degrees = rot_val * i

        # Add a zero as a placeholder
        integrated_signal.append(0)
        integrated_signal_star.append(0)

        # Call the function to get the filtered star and planet spectra
        planet_spectra, star_spectra = filter_spectra(bb_spectra, star_spectra, filt,
                                                        transmission_filter_name, wav_subset, phase_degrees)

        # Get the star and planet plus star signals
        signal, signal_star = get_fp_fs(interp_function, star_radius, planet_spectra, planet_radius, star_spectra)

        # Integrate the total signal and the star signal
        integrated_signal[i] = trapz(signal, x=planet_spectra.wavelength * 1e6)
        integrated_signal_star[i] = trapz(signal_star, x=planet_spectra.wavelength * 1e6)

    # Convert both the signal and the star signal to arrays
    integrated_signal = np.asarray(integrated_signal)
    integrated_signal_star = np.asarray(integrated_signal_star)

    # Divide out the two integrated signals
    fp_fs_ratio = (integrated_signal / integrated_signal_star - 1.0)
    
    # Save the data
    #pd.DataFrame({'Phase': np.arange(0, 360, rot_val), 'Fp_Fs_pmm': fp_fs_ratio * 1e6}
    #                ).to_csv('OUTPUT_DATA/Blackbody_Phase_Curves.txt'.format, sep=' ')
    
    """
    # Plot the data
    phases = np.linspace(0, 345, num_phases) / 360
    ax.plot(phases, fp_fs_ratio * 1e6,
            #linestyle='solid',
            #color=my_colors(k / len(planet_names)),
            linewidth=2,
            label=planet_name)

    # Figure legend stuff
    ax.set_xlim(min(phases), max(phases))
    ax.legend(fontsize=12, loc=(0, 1.03), ncol=2, mode='expand')
    ax.set_xlabel('Orbital Phase')
    ax.set_ylabel(r'F$_p$/F$_s$ (ppm)')
    plt.savefig('../Figures/Blackbody_Phase_Curve.jpg', dpi=200, bbox_inches='tight')
    #plt.clf()
    """
    return None
