"""
Plot a chain from the MCMC run.
Plots values as a function of step in the chain, then makes a corner plot
of the resutls after the burn-in portion of the chain.
"""

import numpy as np
import matplotlib.pyplot as plt
import zmcwrap

def plotChain(fileName, nBurn, par_exclude=None):
    """
    Read a chain file, and cut the first nBurn steps in the emcee run.
    Plot the chain as a function of step in the chain, and make a 
    corner plot of the good portion of the chain (after removing nBurn steps).
    Optionally, don't plot parameters in the list of names from par_exclude.
    """
    
    #Read chain file
    labels, samples, npars, nsteps, nwalkers, fixedParams, wlRanges \
        = zmcwrap.readChain(fileName)

    # Remove parameters that are not to be plotted
    if par_exclude is not None:
        print('excluding from the plots:', par_exclude)
        for par in par_exclude:
            if par in labels:
                exclude_index = labels.index(par)
                labels.remove(par)
                samples = np.delete(samples, exclude_index, axis=1)
                npars -= 1

    # Rename abun_## labels to atomic symbols
    elements=('H' ,'He','Li','Be','B' ,'C' ,'N' ,'O' ,'F' ,'Ne','Na','Mg',
              'Al','Si','P' ,'S' ,'Cl','Ar','K' ,'Ca','Sc','Ti','V' ,'Cr',
              'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
              'Rb','Sr','Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
              'In','Sn','Sb','Te','I' ,'Xe','Cs','Ba','La','Ce','Pr','Nd',
              'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf',
              'Ta','W' ,'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po',
              'At','Rn','Fr','Ra','Ac','Th','Pa','U' )
    for i, label in enumerate(labels):
        if label[:5] == 'abun_':
            atom_num = int(label[5:])
            print('converting ', label, ' to ', elements[atom_num - 1])
            labels[i] = elements[atom_num - 1]
        if (label[:6] == 'contN_' or label[:6] == 'contC_'
            or label[:6] == 'contL_'):
            labels[i] = 'c_'+label[6:]
    
    # chainPlot wants a 3D array with shape (npars,nsteps,nwalkers)
    parsByWalker = samples.T.reshape((npars,nsteps,nwalkers))
    
    # Plot chains
    fig = zmcwrap.chainPlot(parsByWalker, labels)
    plt.show()
    
    # Cut burn-in part of the chain
    if nBurn < 0: nBurn = nsteps//4
    print('trimming the first {:} steps'.format(nBurn))
    samplesT = samples[nwalkers*nBurn:, :]
    # And print the medians and uncertainties
    zmcwrap.zmcwrapMain.printFinalParams(samplesT, labels, fixedParams)
    
    # Make a corner plot
    #import corner #If you want to use the corner package
    #fig = corner.corner(samplesT, labels=labels, bins=20)
    fig = zmcwrap.cornerPlotLight(samplesT, labels, bins=20,
                                  contour_sigmas=[1.0, 2.0],
                                  hist_sigmas=[1.0, 2.0],
                                  figsize=(10, 10))
    fig.savefig('corner.pdf')
    plt.show()
    return

#For running this script as a command line program
if __name__ == "__main__":
    #Optional input files as command line arguments
    import argparse
    parser = argparse.ArgumentParser(description='Plot a chain from an MCMC run.')
    parser.add_argument('file', nargs='?', default='chain.dat',
                        help='Input file with the chain from the MCMC run. Optional, defaults to chain.dat')
    parser.add_argument('-b', '--burn', default=-1, type=int,
                        help='Number of MCMC steps to burn, removing them from the final analysis.')
    parser.add_argument('-ep', '--exclude_par', nargs='*', default=None,
                        help="Parameter(s) to exclude from the corner plot and chain plot. Multiple parameters can be given here, seperated by spaces.")
    args = parser.parse_args()
    fileName = args.file
    nBurn = args.burn
    par_exclude = args.exclude_par
    if isinstance(par_exclude, str): par_exclude = par_exclude.split()

    plotChain(fileName, nBurn, par_exclude)
