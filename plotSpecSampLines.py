"""
Plot synthetic spectra for sampels from an MCMC chain,
along with an observation.
"""

import numpy as np
import matplotlib.pyplot as plt
import specpolFlow as pol
import zmcwrap

def plotSpecSamp(obsFile, synBestFile, synSamplesFile, lineListFile,
                 depthCutoff=0.05, wlMin=None, wlMax=None, figsize=(10, 6)):
    
    specSet = zmcwrap.readSpecSamples(synSamplesFile)
    nspec = specSet.shape[0] - 1
    
    if wlMin is None: wlMin = specSet[0, 0]
    if wlMax is None: wlMax = specSet[0, -1]
    if wlMax < specSet[0, 0] or wlMin > specSet[0, -1]:
        print('WARNING: given wavelength range outside model spectrum range!')
    
    fig, ax = plt.subplots(layout='constrained', figsize=figsize)
    
    # Plot the set of sample synthetic spectra
    for i in range(nspec):
        ax.plot(specSet[0,:], specSet[i+1,:], color='blue', alpha=0.02, lw=2.0)
    
    # Plot the observation
    obsWl, obsSpec, obsErr = zmcwrap.readObservation.getObsSpec(obsFile)
    #plt.errorbar(obsWl, obsSpec, yerr=obsErr, color='0.7', zorder=-2, lw=1.0)
    ax.plot(obsWl, obsSpec, color='k', lw=1.0, zorder=-1)
    
    # Plot the best model
    specBest = np.loadtxt(synBestFile, unpack=True)
    ax.plot(specBest[0,:], specBest[1,:], color='cyan', lw=1.0, zorder=10)
    
    ax.set_xlabel('Wavelength (A)')
    ax.set_ylabel('Flux')

    # Plot line list tick marks and labels
    lineList = pol.read_VALD(lineListFile)
    pol.plot_LineList(lineList, depthCut=depthCutoff,
                      scaleDepths=0.2, rise=0.02, maxLabels=100, ax=ax)
    
    # Set plot view limits
    ax.set_xlim(wlMin, wlMax)
    ymin = np.min(specBest[1, (specBest[0,:] > wlMin) & (specBest[0,:] < wlMax)])
    ymax = 1.04
    ax.set_ylim(ymin - 0.10*(ymax-ymin), ymax + 0.05*(ymax-ymin))

    return fig


#For running this script as a command line program
if __name__ == "__main__":
    #Optional input files as command line arguments
    import argparse
    parser = argparse.ArgumentParser(description='Plot a set of model spectra from an MCMC run, along with an obesrvation.')
    parser.add_argument('-o', '--observation', default='observed.dat',
                        help='Observed spectrum file name')
    parser.add_argument('-b', '--best', default='outSpeci.dat',
                        help='Best fit synthetic spectrum file name')
    parser.add_argument('-s', '--samples', default='outSpecSamples.dat',
                        help='Sampel synthetic spectra file name')
    parser.add_argument('-l', '--lineList', default='data/vlines.dat',
                        help='Line list to plot for line identifications, in a VALD3 long format')
    parser.add_argument('-d', '--depthCutoff', default=0.05, type=float,
                        help='Best fit synthetic spectrum file name')
    parser.add_argument('-w', '--wavelength', default=(None, None),
                        nargs=2, type=float,
                        help='Start and end wavelengths for the range plotted')
    args = parser.parse_args()
    obsFile = args.observation
    synBestFile = args.best
    synSamplesFile = args.samples
    lineListFile = args.lineList
    depthCutoff = args.depthCutoff
    wlMin, wlMax = args.wavelength
    
    fig = plotSpecSamp(obsFile, synBestFile, synSamplesFile, lineListFile, 
                       depthCutoff, wlMin, wlMax)

    fig.savefig('fig-specSet.pdf')
    plt.show()
