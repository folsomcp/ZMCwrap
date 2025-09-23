import numpy as np

from . import readObservation as readObs
from . import zmcwrapMain as mf

def calcSynSpecSamples(nsamples, nthin=0, chainFile='chain.dat',
                       obsFile='observed.dat',
                       specSamplesFile='outSpecSamples.dat'):
    """
    Calculate a set of synthetic spectra, corresponding to samples from the end
    of an MCMC chain, and save them to a file.

    Takes:
    nsamples - The number of samples from the chain to use,
               and the number of synthetic spectra to calculate.
               Starts with the last sample in the chain and works backwards.
    nthin - Optionally, skip every nthin steps in the chain.
            This uses the values for every walker at one step in the chain,
            and then skips nthin for all walkers.
            So the last nwalkers samples will be used, then nthin*nwalkers
            samples will be skipped, then the next nwalkers samples will be
            used, and so on.
    chainFile - The file to read the pre-computed chain from.
    obsFile - The file to read the observation from. This should be consistent
              with what was used to generate chainFile.
              (Only the wavelengths are really used.)
    specSamplesFile - Output file for the set of sample synthetic spectra.  
    """
    
    # Read the chain file, and get both the current parameter values
    # and parameter names.
    freeParID, samples, npars, nsteps, nwalkers, fixedParams, wlRanges\
        = mf.readChain(chainFile)
    
    # Read in the observation, and get the pixels in the range(s) to be fit
    wlObs, specIobs, errObs = readObs.getObsSpecInRange(wlRanges, obsFile)
    
    # Check the number of samples to be used
    nAvailableSamples = samples.shape[0]
    nMaxUsedSamples =  nsamples + ((nsamples-1)//nwalkers)*nwalkers*nthin
    print('using {:} samples from a total of {:}'.format(
        nMaxUsedSamples, nAvailableSamples))
    if nMaxUsedSamples > nAvailableSamples:
        raise ValueError('Requesting more samples than available in the chain! '
                         +'({:} available, requesting {:})'.format(
                             nAvailableSamples, nMaxUsedSamples))
    
    fOut = open(specSamplesFile, 'w')
    for i in range(nsamples):
        # skip every nthin steps,
        # so used nwalkers samples and skip every nwalkers*nthin samples
        iuse = i + (i//nwalkers)*nwalkers*nthin
        freePar = samples[-iuse - 1, :]
        
        # Generate a model spectrum
        lnlike, specWl, specI = mf.lnlike(freePar, freeParID, fixedParams, 
                                          wlObs, specIobs, errObs,
                                          returnSpec=True)
        
        # On the first iteration write wavelengths
        if i == 0:
            fOut.write('{:}\n'.format(specWl.size))
            for j in range(specWl.size):
                fOut.write('{:5.4f}\n'.format(specWl[j]))
            fOut.write('\n')
        
        # Write the intensity spectrum
        fOut.write('{:} {:}\n'.format(i, nAvailableSamples - iuse))
        for j in range(specWl.size):
            fOut.write('{:15.8e}\n'.format(specI[j]))
        fOut.write('\n')
    return


def readSpecSamples(fileName):
    """
    Read in a set of sample model spectra from an MCMC run.
    
    Takes:
    fileName - name the file containing the set of spectra.
    
    Returns:
    specSet - An array containing the model spectra, the first column has
              wavelengths for each pixel, the subsequent columns have the
              intensity values for each model spectrum.
    """
    fin = open(fileName, 'r')
    # First get the dimensions of the grid of spectra to read
    line = fin.readline()
    npix = int(line.split()[0])
    nblock = -1 #don't count the first wavelengths block
    i = 1 #one line already read
    # each block in the file has an index at the start and a blank line at the end
    for line in fin:
        i += 1
        if i == npix + 2:
            nblock += 1
            i = 0
    # Second, read the set of spectra into an array
    # (with the fist column being wavelength)
    specSet = np.zeros((nblock + 1, npix))
    i = 0
    ib = 0
    fin.seek(0)
    for line in fin:
        if i > 0 and i <= npix:
            specSet[ib, i - 1] = float(line.split()[0])
        elif i > npix:
            i = -1
            ib += 1
        i += 1            
    return specSet

