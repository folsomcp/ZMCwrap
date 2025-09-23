"""
Functions for processing observation files.
This can be pretty simple for text format observation files.
"""

import numpy as np

def getObsSpecInRange(wlRanges, obsName='observed.dat'): 
    """
    Read in the observation file for fitting/calculating likelihoods,
    Then return only regions of it in the wavelength ranges specified by 
    wlRanges. wlRanges should be a 2D array (or list of lists),
    with the 1st dimension being different wavelength ranges, 
    and the second dimension being a start and end wavelength.
    """
    _wlRanges = wlRanges
    if isinstance(wlRanges, list) or isinstance(wlRanges, tuple):
        _wlRanges = np.array(wlRanges)
    
    wlObsFull, specIobsFull, errObsFull = getObsSpec(obsName)
    # Sort by wavelength
    # (later interpolation functions may need data sorted by wavelength)
    try:
        # More efficient for spectral points which are mostly already sorted
        isort = np.argsort(wlObsFull, kind='stable')
    except ValueError:
        # But any sort algorithm should be good enough!
        isort = np.argsort(wlObsFull)
    wlObsFull = wlObsFull[isort]
    specIobsFull = specIobsFull[isort]
    errObsFull = errObsFull[isort]
    # Get just the portions of the observation to fit
    indObs = sliceObsWlRange(wlObsFull, _wlRanges)
    wlObs = wlObsFull[indObs]
    specIobs = specIobsFull[indObs]
    errObs = errObsFull[indObs]
    return wlObs, specIobs, errObs

def getObsSpec(obsName='observed.dat'):
    """
    Read in the observation file for fitting/calculating likelihoods.
    """
    wl = []
    specI = []
    err = []
    fIn = open(obsName, 'r')
    line = fIn.readline()
    ncols = len(line.split())
    fIn.seek(0)
    for line in fIn:
        vals = line.split()
        # For 6 column polarimetric LibreESPRIT .s format
        if ncols == 6:
            wl.append(float(vals[0]))
            specI.append(float(vals[1]))
            err.append(float(vals[5]))
        # For 3 column intensity .s format
        elif ncols == 3:
            wl.append(float(vals[0]))
            specI.append(float(vals[1]))
            err.append(float(vals[2]))
        else:
            raise ValueError('Got unexpected number of columns when reading '
                             'observation file {:}, {:}'.format(obsName, ncols))
    fIn.close()
    return np.array(wl), np.array(specI), np.array(err)


def sliceObsWlRange(wlObs, wlRanges):
    """
    For each start and end pair of wavelengths, get the parts of the 
    observation in that range.  wlRanges should be a 2D array, with
    the 1st dimension being different wavelength ranges,
    and the second dimension being a start and end wavelength.
    """
    indUseAll = np.zeros_like(wlObs, dtype=bool)
    for i in range(wlRanges.shape[0]):
        indUseAll += (wlObs >= wlRanges[i,0]) & (wlObs <= wlRanges[i,1])
    #Sanity check for wavelength ranges
    if np.all(np.logical_not(indUseAll)):
        raise ValueError ('ERROR: no observed points in requested wavelength '
                          'range!')
    return indUseAll
