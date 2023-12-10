import numpy as np
from PyAstronomy import pyasl
from scipy.optimize import curve_fit
import matplotlib.pylab as plt
import pandas as pd
import matplotlib.colors as mcolors

# colormap
cm_name = 'romaO'
cm_file = np.loadtxt(f'ScientificColourMaps8/{cm_name}/{cm_name}.txt')
my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file[::-1])

def get_cc(file_dopp, file_rest):
    df_dopp = pd.read_csv(file_dopp, delim_whitespace=True, names=['wav', 'flux', 'reflected'])
    df_rest = pd.read_csv(file_rest, delim_whitespace=True, names=['wav', 'flux', 'reflected'])

    w = np.asarray(list(df_dopp.wav))
    f = np.asarray(list(df_dopp.flux))

    tw = np.asarray(list(df_rest.wav))
    tf = np.asarray(list(df_rest.flux))

    # Carry out the cross-correlation.
    # The RV-range is -10 to +10 km/s in steps of 0.1 km/s.
    # The first and last 200 points of the data are skipped.
    rv, cc = pyasl.crosscorrRV(w, f, tw, tf, -15, 15, 0.1, skipedge=200)
    
    # normalize cc function 
    cc = (cc - min(cc)) / (max(cc) - min(cc))
        
    # Find the index of maximum cross-correlation function
    maxind = np.argmax(cc)

    # fit Guassian to peak region of cc function
    rv_range = rv[maxind - 20: maxind + 20]
    cc_range = cc[maxind - 20: maxind + 20]

    popt, pcov = curve_fit(gaussian, rv_range, cc_range, p0=[1, 0, 1])

    # calculate best fit guassian
    rv_fit = np.linspace(rv[maxind - 20], rv[maxind + 20], 1000)
    cc_fit = gaussian(rv_fit, *popt)

    # index of maximum of Gaussian
    max_gauss = np.argmax(cc_fit)

    return rv, cc, rv_fit, cc_fit, max_gauss

# define Gaussian function
def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0) * (x - x0) / (2 * sigma * sigma))


def plot_cross_correlations(planet_names, num_phases):
    base_file = '../Spectral-Processing/FINISHED_SPECTRA/Spec_{}_{}_phase_{}_inc_0.00.00.0000.00.dat'

    temp = list(np.linspace(0, 360, num_phases, endpoint=False))
    phases = [str(i) for i in temp]
    colors = np.linspace(0, 1, len(phases) + 1)

    doppler_shifts = {}

    for planet_name in planet_names:
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))

        doppler_shifts[planet_name] = []
        for i, phase in enumerate(phases):
            # Get the files based on the phases
            file_dopp = base_file.format(1, planet_name, phase)
            file_rest = base_file.format(0, planet_name, phase)

            # Get the cross correlation
            rv, cc, rv_fit, cc_fit, max_gauss = get_cc(file_dopp, file_rest)

            # plot CCFs
            doppler_shifts[planet_name].append(rv_fit[max_gauss])

            ax.plot(rv, cc + i, lw=2, color=my_colors(colors[i]))
            ax.plot(rv_fit[max_gauss], cc_fit[max_gauss] + i, '.', color=my_colors(colors[i]), ms=15)
            
        sm = plt.cm.ScalarMappable(cmap=my_colors, norm=plt.Normalize(vmin=0, vmax=1))
        sm._A = []
        cbar = fig.colorbar(sm, aspect=30, pad=0.02)
        cbar.set_label('Orbital phase', fontsize=22)

        ax.set_ylabel('Cross Correlation + Offset')
        ax.set_xlabel('Doppler Shift (km/s)')
        ax.set_xlim(-10, 10)
        fig.savefig('../Figures/cc_' + planet_name + '.jpg', dpi=100, bbox_inches='tight')
        plt.clf()


    fig, ax = plt.subplots(1, 1, figsize=(16.18, 10))
    plt.subplots_adjust(wspace=0.06, hspace=0)

    colors = ['black', 'purple', 'green',
              'black', 'purple', 'green']
    
    linestyle_strs = ['solid', 'solid', 'solid',
                      'dashed', 'dashed', 'dashed']
              


    for i, planet_name in enumerate(planet_names):
        final_doppler_shifts = doppler_shifts[planet_name] + [doppler_shifts[planet_name][0]]
        ax.plot(np.linspace(0, 1.0, len(phases) + 1), final_doppler_shifts,label=planet_name)#,
              #  linewidth=1.5, linestyle=linestyle_strs[i], color=colors[i])
    ax.axhline(y=0.0, color='black', linestyle='dotted', linewidth=1.5)

    ax.legend(fontsize=12, loc=(0, 1.03), ncol=3, mode='expand', title_fontsize=18)
    ax.set_ylim(-10, 10)
    ax.set_xlabel('Orbital Phase')
    ax.set_ylabel('Doppler Shift (km/s)')
    ax.set_xlim(0, 1.0)
    fig.savefig('../Figures/Doppler_shifts_' + planet_name + '.jpg', dpi=100, bbox_inches='tight')