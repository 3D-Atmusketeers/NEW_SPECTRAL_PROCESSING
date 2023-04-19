#!/usr/bin/env python
# coding: utf-8

import pandas as pd

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


def filter_spectra(planet_spectra, star_spectra, filt, transmission_filter_name, wav_subset, phase_degrees):
    # Load in the planet spectra and convert to per microns
    #planet_spectra = pd.read_csv(file_path.format(str(phase_degrees)), header=None, delim_whitespace=True, names=['wavelength', 'flux', 'reflected'])
    #planet_spectra.flux = planet_spectra.flux * (3.0e8 / planet_spectra.wavelength ** 2) / 1e6

    if transmission_filter_name == 'MIRI':
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
        filter_5_12 = pd.read_csv('DATA/MIRI_BANDPASS.txt', delim_whitespace=True, skiprows=1, names=['wav', 'transmission'])
        filter_5_12.wav = filter_5_12.wav / 1e6
        f_5_12 = interpolate.interp1d(filter_5_12.wav, filter_5_12.transmission)
        return filter_5_12, f_5_12
    else:
        return None, None


def get_star_spectra(planet_name):
    if ("GJ1214".lower() in planet_name.lower()):
        print("Using a GJ1214 stellar spectrum")
        # Get the initial flux
        stellar_spectrum_1 = pd.read_csv(
            'DATA/GJ1214b_stellar_spectrum.txt',
            names=[
                'wavelength',
                'flux_Watt/m^2/micron'],
            delim_whitespace=True,
            skiprows=1)
        wavelengths_meters = np.asarray(list(stellar_spectrum_1.wavelength * 1e-6)) * u.m
        flux_Watt_m_2_microns = np.asarray(list(stellar_spectrum_1['flux_Watt/m^2/micron'])) * u.W / u.m**2 / u.um
        flux_si = flux_Watt_m_2_microns * np.pi * u.sr
        star_spectra = pd.DataFrame(list(zip(wavelengths_meters.value, flux_si.value)), columns=['wavelength', 'flux'])

        star_radius = 0.204 * 6.957e8

    elif ("peg51".lower() in planet_name.lower()):
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
    elif ("Taub".lower() in planet_name.lower()):
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
    elif ("HD189".lower() in planet_name.lower()):
        print("Using a HD189 stellar spectrum")

        # Get the initial flux
        # This is in nm and ergs/s/cm**2/ster/nm
        stellar_spectrum = pd.read_csv('DATA/HD189733b_stellar_spectrum.txt',
                                       names=['wavelength', 'flux', 'hmm'],
                                       delim_whitespace=True)
        wavelengths_meters = np.asarray(list(stellar_spectrum.wavelength * 1e-9)) * u.m

        # Interestingly, you get 1e-7 from ergs/s to Watts, and 1e4 * 1e3 from cm to m and nm to um
        flux_Watt_m_2_microns = np.asarray(list(stellar_spectrum['flux'])) * u.W / u.m**2 / u.um
        flux_si = flux_Watt_m_2_microns * u.sr * 4.0 * np.pi
        star_spectra = pd.DataFrame(list(zip(wavelengths_meters.value, flux_si.value)), columns=['wavelength', 'flux'])

        star_radius = 0.805 * 6.957e8
    elif ("HD209".lower() in planet_name.lower()):
        print("Using a HD209 stellar spectrum")

        # Get the initial flux
        # This is in nm and ergs/s/cm**2/ster/nm
        stellar_spectrum = pd.read_csv(
            'DATA/HD209458b_stellar_spectrum.txt',
            names=['wavelength', 'flux', 'hmm'],
            delim_whitespace=True)
        wavelengths_meters = np.asarray(list(stellar_spectrum.wavelength * 1e-9)) * u.m

        # Interestingly, you get 1e-7 from ergs/s to Watts, and 1e4 * 1e3 from cm to m and nm to um
        flux_Watt_m_2_microns = np.asarray(list(stellar_spectrum['flux'])) * u.W / u.m**2 / u.um
        flux_si = flux_Watt_m_2_microns * u.sr * 4.0 * np.pi
        star_spectra = pd.DataFrame(list(zip(wavelengths_meters.value, flux_si.value)), columns=['wavelength', 'flux'])

        star_radius = 1.203 * 6.957e8
    else:
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
def plot_planet_spectra_blackbody_comparison_hz(planet_names, black_body_temperatures, num_phases):
    """
    plot the planet spectra, compared to a blackbody
    planet names: list of strings
    black_body_temperatures: list of floats
    """

    colors_bb = pl.cm.Greys(np.linspace(0, 1, len(black_body_temperatures) + 2))

    cm_name = 'batlow'
    cm_file = np.loadtxt(f'ScientificColourMaps7/{cm_name}/{cm_name}.txt')
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
            if planet_name == 'Peg51b':
                INC_STRING = '0.17450'
            elif planet_name == 'Taub':
                INC_STRING = '0.7850'
            else:
                INC_STRING = '0.00'

            file_path = '../Spectral-Processing/FINISHED_SPECTRA/Spec_0_' + planet_name + '_phase_{}_inc_' + INC_STRING + '.00.0000.00.dat'

            # Load in the planet spectra
            planet_spectra = pd.read_csv(file_path.format(str(i * rot_val)), header=None,
                                         delim_whitespace=True, names=['wavelength', 'flux', 'reflected'])

            # Reset the index
            planet_spectra = planet_spectra.reset_index(drop=True)

            ax.plot(planet_spectra.wavelength * 1e6,
                    planet_spectra.flux,
                    alpha=1.0,
                    color=my_colors(i / num_phases),
                    linewidth=0.5,
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
        ax.set_ylim(1e-12, 1e-8)
        ax.set_yscale('log')
        ax.set_xlim(5, 20)
        ax.legend(fontsize=12, loc=(0, 1.05), ncol=2, mode='expand', title_fontsize=16)
        ax.set_xlabel(r'Wavelength ($\mu$m)')
        ax.set_ylabel(r'Flux (W/m$^2$Hz)')  # (W m$^{-2}$)
        plt.savefig('../Figures/planet_spectra_blackbody_comparison_hz_{}.jpg'.format(planet_name), dpi=200, bbox_inches='tight')
        plt.clf()
    return None



def plot_planet_spectra_blackbody_comparison_microns(planet_names, black_body_temperatures, num_phases):
    """
    plot the planet spectra, compared to a blackbody
    planet names: list of strings
    black_body_temperatures: list of floats
    """

    colors_bb = pl.cm.Greys(np.linspace(0, 1, len(black_body_temperatures) + 2))

    cm_name = 'batlow'
    cm_file = np.loadtxt(f'ScientificColourMaps7/{cm_name}/{cm_name}.txt')
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
            if planet_name == 'Peg51b':
                INC_STRING = '0.17450'
            elif planet_name == 'Taub':
                INC_STRING = '0.7850'
            else:
                INC_STRING = '0.00'

            # Get the file path
            file_path = '../Spectral-Processing/FINISHED_SPECTRA/Spec_0_' + planet_name + '_phase_{}_inc_' + INC_STRING + '.00.0000.00.dat'

            # Load in the planet spectra
            planet_spectra = pd.read_csv(file_path.format(str(i * rot_val)), header=None,
                                         delim_whitespace=True, names=['wavelength', 'flux', 'reflected'])

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
        #ax.set_ylim(1e-12, 1e-8)
        ax.set_ylim(0, 2000)
        ax.set_xlim(5, 20)
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
        Teffs = [3000, 4000, 5000, 6000]
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
                linewidth=1.5)

        plt.ylim(1e4, 3e8)
        plt.xlim(0.1, np.max(star_spectra.wavelength * 1e6))
        plt.yscale('log')
        plt.xlabel('Wavelengths (microns)')
        plt.ylabel(r'Flux (Watts/m$^2$/micron)')
        plt.legend(loc='upper right')
        plt.savefig('../Figures/star_test_{}.jpg'.format(planet_name), dpi=200, bbox_inches='tight')
        plt.clf()
    return None


def plot_filters(planet_names):
    filter_5_12, f_5_12 = get_filter(which_filter='MIRI')

    

    # Figure aesthetics
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
    plt.subplots_adjust(hspace=0.05, wspace=0.25)

    # Get the file path
    if planet_names[0] == 'Peg51b':
        INC_STRING = '0.17450'
    elif planet_names[0] == 'Taub':
        INC_STRING = '0.7850'
    else:
        INC_STRING = '0.00'

    file_path = '../Spectral-Processing/FINISHED_SPECTRA/Spec_0_' + planet_names[0] + '_phase_{}_inc_' + INC_STRING + '.00.0000.00.dat'
    spectra = pd.read_csv(
        file_path.format(
            str(0.0)),
        header=None,
        delim_whitespace=True,
        names=[
            'wavelength',
            'flux',
            'reflected'])

    # Get only the spectra in the filter ranges we care about
    #spectra_3_6  = spectra[(spectra['wavelength'] > np.amin(filter_3_6.wav)) & (spectra['wavelength'] < np.amax(filter_3_6.wav))]
    #spectra_4_5  = spectra[(spectra['wavelength'] > np.amin(filter_4_5.wav)) & (spectra['wavelength'] < np.amax(filter_4_5.wav))]
    spectra_5_12 = spectra[(spectra['wavelength'] > np.amin(filter_5_12.wav)) & (spectra['wavelength'] < np.amax(filter_5_12.wav))]


    # Plot the interpolated transmission stuff and the original ones
    # ax.plot(filter_3_6.wav, filter_3_6.transmission, color='#7e1e9c', linewidth=3, alpha=0.2, label="")
    # ax.plot(filter_4_5.wav, filter_4_5.transmission, color='#15b01a', linewidth=3, alpha=0.2, label="Original")
    ax.plot(filter_5_12.wav * 1e6, filter_5_12.transmission, color='red', linewidth=3, alpha=0.9, label="MIRI")

    ax.plot(spectra_5_12.wavelength * 1e6, f_5_12(spectra_5_12.wavelength),
            color='black',
            linewidth=3,
            linestyle='dashed',
            label="Interpolated")

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
            if planet_name == 'Peg51b':
                INC_STRING = '0.17450'
            elif planet_name == 'Taub':
                INC_STRING = '0.7850'
            else:
                INC_STRING = '0.00'

            file_path = '../Spectral-Processing/FINISHED_SPECTRA/Spec_0_' + planet_name + '_phase_{}_inc_' + INC_STRING + '.00.0000.00.dat'

            # Load in the planet spectra
            planet_spectra = pd.read_csv(file_path.format(str(i * rot_val)), header=None,
                                         delim_whitespace=True, names=['wavelength', 'flux', 'reflected'])
            planet_spectra.flux = planet_spectra.flux * (3.0e8 / planet_spectra.wavelength ** 2) / 1e6

            # Reset the index
            planet_spectra = planet_spectra.reset_index(drop=True)


            ax.plot(planet_spectra.wavelength * 1e6,
                    planet_spectra.flux,
                    #color=my_colors(i / num_phases),
                    linewidth=1,
                    label=str(np.round(rot_val * i / 360., 3)))




        # Do somet figure stuff
        #ax.set_xlim(min(planet_spectra.wavelength * 1e6), max(planet_spectra.wavelength * 1e6))
        ax.set_xlim(0.5, 12)
        ax.set_ylim(1e2,1e6)
        ax.legend(fontsize=12, loc=(0, 1.03), ncol=5, mode='expand', title='Orbital Phase', title_fontsize=16)
        ax.set_xlabel(r'Wavelength ($\mu$m)')
        ax.set_ylabel(r'Flux (W/m$^2$/micron)')
        ax.set_yscale('log')
        plt.savefig('../Figures/Planet_Simple_Spectra_{}.jpg'.format(planet_name), dpi=200, bbox_inches='tight')
        plt.clf()

    return None


def plot_fp_fs_phase_curves(planet_names, planet_name_char_len, planet_radii, num_phases,
                            transmission_filter_name, wav_subset):
    # Figure aesthetics
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6), sharex=True, sharey=False)
    plt.subplots_adjust(hspace=0.05, wspace=0.25)

    cm_name = 'batlow'
    cm_file = np.loadtxt(f'ScientificColourMaps7/{cm_name}/{cm_name}.txt')
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

        for i in range(num_phases):
            # For reading in the file names
            rot_val = 360. / num_phases
            phase_degrees = rot_val * i

            # Add a zero as a placeholder
            integrated_signal.append(0)
            integrated_signal_star.append(0)

            # Load in the planet spectra and convert to per microns

            # Get the file path
            if planet_name == 'Peg51b':
                INC_STRING = '0.17450'
            elif planet_name == 'Taub':
                INC_STRING = '0.7850'
            else:
                INC_STRING = '0.00'
            file_path = '../Spectral-Processing/FINISHED_SPECTRA/Spec_0_' + planet_name + '_phase_{}_inc_' + INC_STRING + '.00.0000.00.dat'

            # Call the function to get the filtered star and planet spectra
            planet_spectra, star_spectra = filter_spectra(planet_spectra, star_spectra, filt,
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
        pd.DataFrame({'Phase': np.arange(0, 360, rot_val), 'Fp_Fs_pmm': fp_fs_ratio * 1e6}
                     ).to_csv('OUTPUT_DATA/Fp_Fs_{}_Phase_Curves.txt'.format(planet_name), sep=' ')
        
        print(planet_name)

        if '30met' in planet_name.lower():
            linestyle_str='dashed'
        else:
            linestyle_str='solid'

        # Plot the data
        phases = np.linspace(0, 345, num_phases) / 360
        ax.plot(phases, fp_fs_ratio * 1e6,
                #linestyle='solid',
                #color=my_colors(k / len(planet_names)),
                linestyle=linestyle_str,
                linewidth=2,
                label=planet_name)

    # Figure legend stuff
    ax.set_xlim(min(phases), max(phases))
    ax.legend(fontsize=12, loc=(0, 1.03), ncol=2, mode='expand')
    ax.set_xlabel('Orbital Phase')
    ax.set_ylabel(r'F$_p$/F$_s$ (ppm)')
    plt.savefig('../Figures/Fp_Fs_Phase_Curves.jpg', dpi=200, bbox_inches='tight')
    plt.clf()
    return None


def plot_fp_phase_curves(planet_names, planet_name_char_len, num_phases,
                         transmission_filter_name, wav_subset):
    # Figure aesthetics
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6), sharex=True, sharey=False)
    plt.subplots_adjust(hspace=0.05, wspace=0.25)

    cm_name = 'batlow'
    cm_file = np.loadtxt(f'ScientificColourMaps7/{cm_name}/{cm_name}.txt')
    cm_file = np.roll(cm_file, 140, axis=0)
    my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file)

    # I currently only have MIRI coded in
    filt, interp_function = get_filter(which_filter=transmission_filter_name)

    for k, planet_name in enumerate(planet_names):
        # Set the arrays for the integrated signals
        integrated_planet_spectra = []

        # Get the star spectra
        star_spectra, star_radius = get_star_spectra(planet_name)

        for i in range(num_phases):
            # For reading in the file names
            rot_val = 360. / num_phases
            phase_degrees = rot_val * i

            # Add a zero as a placeholder
            integrated_planet_spectra.append(0)

            # Load in the planet spectra and convert to per microns
            # Get the file path
            if planet_name == 'Peg51b':
                INC_STRING = '0.17450'
            elif planet_name == 'Taub':
                INC_STRING = '0.7850'
            else:
                INC_STRING = '0.00'
            file_path = '../Spectral-Processing/FINISHED_SPECTRA/Spec_0_' + planet_name + '_phase_{}_inc_' + INC_STRING + '.00.0000.00.dat'

            # Call the function to get the filtered star and planet spectra
            planet_spectra, star_spectra = filter_spectra(planet_spectra,star_spectra, filt, transmission_filter_name, wav_subset, phase_degrees)

            if interp_function is None:
                planet_spectra.flux = planet_spectra.flux
            else:
                planet_spectra.flux = planet_spectra.flux * interp_function(planet_spectra.wavelength)

            # Integrate the total planet planet flux
            integrated_planet_spectra[i] = trapz(planet_spectra.flux, x=planet_spectra.wavelength * 1e6)

        # Convert both the signal and the star signal to arrays
        integrated_planet_spectra = np.asarray(integrated_planet_spectra)

        # Divide out the two integrated signals
        fp = integrated_planet_spectra

        # Save the data
        pd.DataFrame({'Phase': np.arange(0, 360, rot_val), 'fp W/m2/micron': fp}
                     ).to_csv('OUTPUT_DATA/Fp_{}_Phase_Curves.txt'.format(planet_name), sep=' ')

        # Plot the data
        phases = np.linspace(0, 345, num_phases) / 360
        ax.plot(phases, fp,
                linestyle='solid',
                #color=my_colors(k / len(planet_names)),
                linewidth=2.0,
                label=planet_name)

    # Figure legend stuff
    ax.set_xlim(min(phases), max(phases))
    ax.legend(fontsize=12, loc=(0, 1.03), ncol=2, mode='expand')
    ax.set_xlabel('Orbital Phase')
    ax.set_ylabel(r'Planet Flux (W/m$^2$/micron)')
    plt.savefig('../Figures/Fp_Phase_Curves.jpg', dpi=200, bbox_inches='tight')
    plt.clf()
    return None


def plot_fp_fs_spectra(planet_names, planet_radii, num_phases, transmission_filter_name, wav_subset, resolution):
    cm_name = 'batlow'
    cm_file = np.loadtxt(f'ScientificColourMaps7/{cm_name}/{cm_name}.txt')
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

        for i in range(num_phases):
            rot_val = 360. / num_phases
            phase_degrees = rot_val * i

            # Get the file path
            if planet_name == 'Peg51b':
                INC_STRING = '0.17450'
            elif planet_name == 'Taub':
                INC_STRING = '0.7850'
            else:
                INC_STRING = '0.00'

            file_path = '../Spectral-Processing/FINISHED_SPECTRA/Spec_0_' + planet_name + '_phase_{}_inc_' + INC_STRING + '.00.0000.00.dat'

            # Call the function to get the filtered star and planet spectra
            planet_spectra, star_spectra = filter_spectra(planet_spectra, star_spectra, filt, transmission_filter_name, wav_subset, phase_degrees)

            # Get the planet to star flux ratio
            signal, signal_star = get_fp_fs(interp_function, star_radius, planet_spectra, planet_radius, star_spectra)

            # Divide out the two integrated signals
            fp_fs_ratio = (signal / signal_star - 1.0)

            if resolution > 0:
                fp_fs_ratio = np.asarray(reduceSpectralResolution(list(planet_spectra.wavelength),list(fp_fs_ratio),
                                                                  R_low=resolution, R_high=None, lambda_mid=None, n=4))

            ax.plot(planet_spectra.wavelength * 1e6,
                    fp_fs_ratio * 1e6,
                    color=my_colors(i / num_phases),
                    linewidth=1.5,
                    label=str(np.round(rot_val * i / 360., 3)))

            pd.DataFrame({'Wavelength (microns)': planet_spectra.wavelength * 1e6, 'Fp_Fs_pmm': planet_spectra.flux}
                         ).to_csv('OUTPUT_DATA/Fp_Fs_Spectra_{}_Spectra_{}.txt'.format(str(i * rot_val), planet_name), sep=' ')

        # Figure legend
        ax.set_xlim(min(planet_spectra.wavelength * 1e6), max(planet_spectra.wavelength * 1e6))
        ax.legend(fontsize=12, loc=(0, 1.03), ncol=5, mode='expand', title='Orbital Phase', title_fontsize=18)
        ax.set_xlabel(r'Wavelength ($\mu$m)')
        ax.set_ylabel(r'F$_p$/F$_s$ (ppm)')  # (W m$^{-2}$)
        plt.savefig('../Figures/Fp_Fs_Spectra_{}.jpg'.format(planet_name), dpi=200, bbox_inches='tight')

    return None


def plot_fp_spectra(planet_names, planet_radii, num_phases, transmission_filter_name, wav_subset, resolution):
    cm_name = 'batlow'
    cm_file = np.loadtxt(f'ScientificColourMaps7/{cm_name}/{cm_name}.txt')
    cm_file = np.roll(cm_file, 140, axis=0)
    my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file)

    # I currently only have MIRI coded in
    filt, interp_function = get_filter(which_filter=transmission_filter_name)

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
            if planet_name == 'Peg51b':
                INC_STRING = '0.17450'
            elif planet_name == 'Taub':
                INC_STRING = '0.7850'
            else:
                INC_STRING = '0.00'
            file_path = '../Spectral-Processing/FINISHED_SPECTRA/Spec_0_' + planet_name + '_phase_{}_inc_' + INC_STRING + '.00.0000.00.dat'

            # Load in the planet spectra
            planet_spectra = pd.read_csv(file_path.format(str(i * rot_val)), header=None,
                                         delim_whitespace=True, names=['wavelength', 'flux', 'reflected'])
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

            ax.plot(planet_spectra.wavelength * 1e6,
                    flux,
                    color=my_colors(i / num_phases),
                    linewidth=1.5,
                    label=str(np.round(rot_val * i / 360., 3)))

            pd.DataFrame({'Wavelength (microns)': planet_spectra.wavelength * 1e6, 'Planet_Flux_micron': flux}
                         ).to_csv('OUTPUT_DATA/Fp_Spectra_{}_{}.txt'.format(str(i * rot_val), planet_name), sep=' ')

        ax.set_xlim(min(planet_spectra.wavelength * 1e6), max(planet_spectra.wavelength * 1e6))
        ax.legend(fontsize=12, loc=(0, 1.03), ncol=4, mode='expand', title='Orbital Phase', title_fontsize=18)
        ax.set_xlabel(r'Wavelength ($\mu$m)')
        ax.set_ylabel(r'F$_p$ (W/m$^2$/micron)')  # (W m$^{-2}$)
        plt.savefig('../Figures/Fp_Spectra_{}.jpg'.format(planet_name), dpi=200, bbox_inches='tight')

    return None


def plot_blackbody_phase_curve(planet_name, planet_radii,num_phases,transmission_filter_name, wav_subset,resolution,temp):
    cm_name = 'batlow'
    cm_file = np.loadtxt(f'ScientificColourMaps7/{cm_name}/{cm_name}.txt')
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
    
    print(fp_fs_ratio)
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