"""
Plot synthetic spectra for sampels from an MCMC chain,
along with an observation.
"""

import numpy as np
import matplotlib.pyplot as plt
import zmcwrap

def plotSpecSamp(obsFile, synBestFile, synSamplesFile):
    
    specSet = zmcwrap.readSpecSamples(synSamplesFile)
    nspec = specSet.shape[0] - 1
    
    fig, ax = plt.subplots(layout='constrained')
    
    # Plot the set of sample synthetic spectra
    for i in range(nspec):
        ax.plot(specSet[0,:], specSet[i+1,:], color='blue', alpha=0.02, lw=2.0)
    #Autoscale the axes on this dataset, but don't autoscale on subsequent datasets
    ax.autoscale_view()
    ax.autoscale(False)
    
    # Plot the observation
    obsWl, obsSpec, obsErr = zmcwrap.readObservation.getObsSpec(obsFile)
    #plt.errorbar(obsWl, obsSpec, yerr=obsErr, color='0.7', zorder=-2, lw=1.0)
    ax.plot(obsWl, obsSpec, color='k', lw=1.0, zorder=-1)
    
    # Plot the best model
    specBest = np.loadtxt(synBestFile, unpack=True)
    ax.plot(specBest[0,:], specBest[1,:], color='r', lw=1.0, zorder=10)
    
    ax.set_xlabel('Wavelength')
    ax.set_ylabel('Flux')
    
    plt.show()
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
    args = parser.parse_args()
    obsFile = args.observation
    synBestFile = args.best
    synSamplesFile = args.samples

    plotSpecSamp(obsFile, synBestFile, synSamplesFile)
