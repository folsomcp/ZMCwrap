import numpy as np
import scipy.special as sf
import emcee

from . import wrapper_zeeman as sspec
from . import readObservation as readObs

#Data likelyhood P(y|x,sigma,params)
#Calls the main modelling program
def lnlike(freePar, freeParID, fixedParams, wlObs, specIobs, errObs,
           returnSpec=False, verbose=False):
    """
    Calculate the log likelyhood ln(P(y|x,sigma,params))
    This runs the main spectrum synthesis code (usually Zeeman),
    and uses an observation to calculate the likelyhood.
    The log likelyhood (assuming uncorrelated Gaussian uncertainties)
    is essentially chi^2 (apart from a normalizing constant).
    
    Takes:
    freePar - list of free parameter values
    freeParID - list of names for those free parameters
    fixedParams - dictionary of fixed parameter names and values
    wlObs - an array of observed spectrum wavelengths
    specIobs - an array of observed spectrum intensities
    errObs - an array of observed spectrum intensity uncertainties
    returnSpec - a flag for returning the model spectrum,
                 as two more arrays of wavelength and flux.
    verbose - a flag for printing extra diagnostic information
              and saving model spectrum

    Returns:
    log likelyhood value

    Takes the free parameters, keywords for the free parameters,
    as well as a dict for the fixed parameters, all passed in to Zeeman.
    Parameters not specified as free or fixed should be read from the standard
    Zeeman input files.
    """
    
    # Calculate the model spectrum for this set of parameters
    wlSyn, specIsyn = sspec.runSpecSynth(freePar, freeParID, fixedParams, verbose)

    # Interpolate the model onto the observation
    synIinterp = np.interp(wlObs, wlSyn, specIsyn)

    errLogF = 0.0
    # Include an extra continuum normalization if necessary
    includeContPoly = False
    contPoly = np.zeros_like(synIinterp)
    wlRef =  0.5*(wlObs[0]+wlObs[-1])
    # Normaize wavelength to be in the range [-1, 1]
    wlNorm = (wlObs - wlRef)/(wlObs[-1] - wlRef)
    for i, key in enumerate(freeParID):
        if key[:5] == 'contN':
            nContPol = int(key[6:])
            contPoly += freePar[i]*(wlObs-wlRef)**nContPol
            includeContPoly = True
        if key[:5] == 'contC':
            #Chebyshev polynomial of the first kind
            nContPol = int(key[6:])
            contPoly += freePar[i]*sf.eval_chebyt(nContPol, wlNorm)
            includeContPoly = True
        if key[:5] == 'contL':
            #Legendre polynomial
            nContPol = int(key[6:])
            contPoly += freePar[i]*sf.eval_legendre(nContPol, wlNorm)
            includeContPoly = True
        if key == 'errLogF':
            # Extra error term
            errLogF = freePar[i]
    for key, value in fixedParams.items():
        if key[:5] == 'contN':
            nContPol = int(key[6:])
            contPoly += value*(wlObs-wlRef)**nContPol
            includeContPoly = True
        if key[:5] == 'contC':
            #Chebyshev polynomial of the first kind
            nContPol = int(key[6:])
            contPoly += value*sf.eval_chebyt(nContPol, wlNorm)
            includeContPoly = True
        if key[:5] == 'contL':
            #Legendre polynomial
            nContPol = int(key[6:])
            contPoly += value*sf.eval_legendre(nContPol, wlNorm)
            includeContPoly = True
    if includeContPoly:
        synIinterp = synIinterp * contPoly

    # Save the model spectrum (usualy a final model)
    if verbose:
        fOuti = open('outSpeci.dat','w')
        for i in range(wlObs.size):
            fOuti.write('{:5.4f} {:13.10f} {:13.10f}\n'.format(
                wlObs[i], synIinterp[i], contPoly[i]))
        fOuti.close()

    #Calculate the log likelyhood (and chi squared)
    if verbose:
        chi2 = np.sum(((specIobs - synIinterp)/errObs)**2)
        print('before scaling uncertainties: chi2', chi2, 'reduced chi2',
              chi2/(wlObs.size-len(freePar)))

    err = errObs*np.exp(errLogF)
    lnlike = -0.5*np.sum( ((specIobs - synIinterp)/err)**2 + np.log(2.*np.pi*err**2) )
    
    if verbose:
        chi2 = np.sum(((specIobs - synIinterp)/err)**2)
        print('with scaled uncertinties: chi2', chi2, 'reduced chi2', 
              chi2/(wlObs.size-len(freePar)), 'lnlike', lnlike)
    
    if returnSpec:
        return lnlike, wlObs, synIinterp
    return lnlike


#Prior probability P(params)
def lnprior(freePar, freeParID):
    """
    Calculate the priors for the parameters ln(P(params))
    
    Takes:
    freePar - a list of parameter values (params)
    freeParID - parameter keyword names (freeParID)
    
    Returns:
    log prior probability (often 0) or -inf for invalid values.
    """
    if(len(freePar) != len(freeParID)):
        raise ValueError ('missmatch between freePar and freeParID lengths')

    ffsum = 0.0
    #Require parameters to be in valid ranges
    for i in range(len(freePar)):
        if freeParID[i] == 'vsini':
            if freePar[i] < 0.0:
                return -np.inf
        if freeParID[i] == 'Vmic':
            if freePar[i] < 0.0:
                return -np.inf
        if freeParID[i] == 'Vmac':
            if freePar[i] < 0.0:
                return -np.inf
        if freeParID[i] == 'Teff':
            if freePar[i] > 30000.0 or freePar[i] < 3500.0:
                return -np.inf
        if freeParID[i] == 'logg':
            if freePar[i] > 5.0 or freePar[i] < 2.5:
                return -np.inf
        if freeParID[i] == 'contFlx':
            if freePar[i] < 0.0:
                return -np.inf
        
        if freeParID[i][:6] == 'FFmono':
            if freePar[i] < 0.0:
                return -np.inf
            #Find some way to enforce filling factors summing to <= 1.0?
            #It could suffice to just reject those models.
            ffsum += freePar[i]
            if ffsum > 1.0:
                return -np.inf
        if freeParID[i][:5] == 'FFdip':
            if freePar[i] < 0.0:
                return -np.inf
            #Reject models with a total filling factor > 1
            ffsum += freePar[i]
            if ffsum > 1.0:
                return -np.inf

        if freeParID[i][:4] == 'abun':
            if freePar[i] > -0.5 or freePar[i] < -13.0:
                return -np.inf
        
    return 0.0


#Full posterior probability P(params|x,y,sigma) ~ P(params)*P(y|x,sigma,params)
def lnprob(freePar, freeParID, fixedParams, wlObs, specIobs, errObs):
    """
    Calculate the full posterior probability P(params|x,y,sigma)
    as ln(P(params|x,y,sigma)) = ln(P(params)) + ln(P(y|x,sigma,params)).
    This relies on the lnprior() and lnlike() functions.
    """
    lnpri = lnprior(freePar, freeParID)
    #If these parameters are bad reject them before running the model.
    if not np.isfinite(lnpri):
        return -np.inf
    lnlik = lnlike(freePar, freeParID, fixedParams, wlObs, specIobs, errObs)
    #If these parameters are off the grid, reject the point too.
    if not np.isfinite(lnlik):
        return -np.inf
    return lnpri + lnlik


def parseInputParams(params, epsilons, fixedParams):
    """
    Conver the input dicts of parameter keywords and values into a list of 
    parameters for emcee and a list of parameter keywords/names for other 
    routines later.  Do this for the free parameter values ('params'), and 
    the initial dispersions of free parameters ('epsions').
    """

    #Check that the parameter keywords are know names
    knownNames = ['Vr', 'vsini', 'Vmic', 'Vmac', 'Teff', 'logg', 'metal',
                  'contFlx', 'Bmono', 'FFmono', 'Bdip', 'FFdip', 'abun_',
                  'contN_', 'contC_', 'contL_', 'errLogF']
    for key in list(params)+list(fixedParams):
        if key[:5] == 'abun_':
            try:
                atom = int(key[5:])
            except:
                raise ValueError('Format error for keyword {:}'.format(key))
        elif key[:6] == 'contN_' or key[:6] == 'contC_' or key[:6] == 'contL_':
            try:
                coeff = int(key[6:])
            except:
                raise ValueError('Format error for keyword {:}'.format(key))
        elif key not in knownNames:
            raise ValueError('Unknown keyword name {:}'.format(key))
    for key in list(params):
        if key in list(fixedParams):
            raise ValueError('Keyword in both free and fixed params '
                             '{:}'.format(key))
    
    #Clean up an possible messiness in the fixedParams dict
    #(mostly things that could be lists or single values: make them all lists)
    print('Fixed parameter values:')
    for key, value in fixedParams.items():
        if key in ('Bmono', 'FFmono', 'Bdip', 'FFdip'):
            if isinstance(value, (int, float)):
                fixedParams[key] = [value]
        print(' ', key, value)
    
    #Build lists of the free parameter values, free parameter keyword IDs,
    #and dispersion for the initial free parameter distribution.
    freePar = []
    freeParID = []
    freeParEps = []
    print('Using free fitting parameters with initial values:')
    for key, value in params.items():
        if value != None:
            if key in epsilons:
                #If the value exists and there is a corresponding epsilon
                print(' ', key, value, '+/-', epsilons[key])
                if isinstance(value, (int, float, complex)):
                    freePar += [value]
                    freeParID += [key]
                    freeParEps += [epsilons[key]]
                elif isinstance(value, (list, tuple)):
                    for val in value:
                        freePar += [val]
                        freeParID += [key]
                    for val in epsilons[key]:
                        freeParEps += [val]
                else:
                    raise ValueError ('ERROR unsupported data type in dict of '
                                'free parameters {:} {:}'.format(key, value))
            else:
                #If there is no correspoinding epsilon
                raise ValueError ('ERROR missing epsilon for free value: '
                                  '{:} {:} '.format(key, value))
            
            # Make sure there are the right number of magnetic field strengths
            # if using filling factors
            if key == 'FFmono':
                if ('Bmono' in params):
                    if len(value) != len(params['Bmono']):
                        raise ValueError ('ERROR: Inconsistent number of Bmono '
                                          'and FFmono values')
                elif ('Bmono' in fixedParams):
                    if len(value) != len(fixedParams['Bmono']):
                        raise ValueError ('ERROR: Inconsistent number of Bmono '
                                          'and FFmono values')
                else:
                    raise ValueError ('ERROR: Missing list of Bmono')
            if key == 'FFdip':
                if ('Bdip' in params):
                    if len(value) != len(params['Bdip']):
                        raise ValueError ('ERROR: Inconsistent number of Bdip '
                                          'and FFdip values')
                elif ('Bdip' in fixedParams):
                    if len(value) != len(fixedParams['Bdip']):
                        raise ValueError ('ERROR: Inconsistent number of Bdip '
                                          'and FFdip values')
                else:
                    raise ValueError ('ERROR: Missing list of Bdip')

    if len(freePar) != len(freeParEps):
        raise ValueError ('mismatch in number of free parameter initial values'\
            +' and epsilons: {:} != {:}'.format(len(freePar), len(freeParEps)))
    return freePar, freeParID, freeParEps


def writeChainFileHeader(freeParID, fixedParams=None, wlRanges=None,
                         fileName="chain.dat"):
    """
    Write out the starting header of the output chain file.
    Takes the set of free parameter names, and optionally the chain file name
    """
    foutChain = open(fileName, "w")
    foutChain.write('#'+(' '.join(freeParID))+'\n')
    if fixedParams is not None:
        outStr = '#'
        for key, value in fixedParams.items():
            outStr += ' {:}: {:},'.format(key, value)
        if len(outStr) > 1: outStr = outStr[:-1]
        foutChain.write(outStr + '\n')
    if wlRanges is not None:
        outStr = '#wl'
        for wlran in wlRanges:
            outStr += ' {:} {:},'.format(wlran[0], wlran[1])
        foutChain.write(outStr[:-1] + '\n')
    foutChain.close()
    return

def writeChainIteration(position, lnprobability=0.0,
                        fileName="chain.dat", verbose=False):
    foutChain = open(fileName, "a")
    for k in range(position.shape[0]): #loop over walkers at this step
        str_values = ""
        for value in position[k]:
            str_values += "{:16e} ".format(value)
        if verbose == True:
            foutChain.write("{:4d} {:s} {:}\n".format(k, str_values,
                                                      lnprobability[k]))
        else:
            foutChain.write("{:4d} {:s}\n".format(k, str_values))
    
    foutChain.close()
    return


def readChain(fileName):
    """
    Read in a text chain file.
    
    Takes:
    fileName - the name of a chain file (may include a path to the file)
    
    Returns:
    labels - list of text strings for the parameter names
    samples - array of parameter values through the chain
    npars - number of model parameters to estimate values
    nsteps - number of steps in the chain
    nwalkers - number of walkers in the chain
    """
    #Read the header line of fitting parameters
    inFile = open(fileName, 'r')
    header = inFile.readline()
    labels = header.strip('# ').split()
    #Read the optional line of fixed parameters
    fixedParams = None
    header = inFile.readline().strip()
    if header[0] == '#':
        fixedParams = {}
        header = header.strip('# ')
        if len(header) > 0:
            parts = header.split(',')
            for part in parts:
                key, val = part.split(':')
                fixedParams[key.strip()] = float(val)
    #Read the optional line of wavelength ranges
    wlRanges = None
    header = inFile.readline().strip()
    if header[0:3] == '#wl':
        wlRanges = []
        parts = header.strip('#wl ').split(',')
        for part in parts:
            start, stop = part.split()
            wlRanges += [[float(start), float(stop)]]
    inFile.close()
    
    #Read the main data
    chainFile = np.loadtxt(fileName, unpack=True)
    nColums = chainFile.shape[0]
    npars = nColums-1
    
    walkerID = chainFile[0,:] 
    samples =  chainFile[1:,:].T
    
    nwalkers = int(np.max(walkerID))+1
    
    walkerID = walkerID.reshape((-1, nwalkers))
    nsteps = walkerID.shape[0]

    if fixedParams is None or wlRanges is None:
        return labels, samples, npars, nsteps, nwalkers
    return labels, samples, npars, nsteps, nwalkers, fixedParams, wlRanges


def printFinalParams(samples, freeParID, fixedParams, 
                     lnprobability=0.0, sampler=None, verbose=False):
    """
    Print some information about the best parameters from the chain,
    (and optionally final step in the chain if verbose=True).
    """
    npars = samples.shape[-1]
    # Print the autocorrelation time and acceptance fraction
    # Using tol=0 means that we'll always get an estimate,
    # even if it isn't trustworthy.
    if sampler is not None:
        print('acceptance fraction {:}'.format(
            np.average(sampler.acceptance_fraction)))
        print('autocorrelation time estimate {:}'.format(
            sampler.get_autocorr_time(tol=0) ))
    
    #Optionally, print output about the final state
    if verbose and sampler is not None:
        print('Final step in the chain, averaged over walkers')
        position = sampler.chain[:, -1, :]
        for i in range(npars):
            print('ending par ({:}) avg {:} std {:} lnP {:}'.format(
                freeParID[i], np.average(position[:,i]), 
                np.std(position[:,i]), np.average(lnprobability) ))

    #Print central and 1sigma values from percentiles
    indFF = 0
    indFF2 = 0
    print('Final best parameters, and 1sigma uncertainties, after burn-in')
    for i in range(npars):
        par_low, par_mid, par_hig, par_m3sig, par_p3sig = np.percentile(samples[:,i],
                                            [15.8655, 50.0, 84.1345, 0.13499, 99.86501])
        if freeParID[i] == 'FFmono':
            pparName = '{:}_{:}'.format(freeParID[i], fixedParams['Bmono'][indFF])
            indFF += 1
        elif freeParID[i] == 'FFdip':
            pparName = '{:}_{:}'.format(freeParID[i], fixedParams['Bdip'][indFF2])
            indFF2 += 1
        else:
            pparName = freeParID[i]
        print('{:7}  {:9.6} {:+10.6} {:+10.6} (3sigma range {:8.6} {:8.6} )'.format(
            pparName, par_mid, par_hig-par_mid, par_low-par_mid, par_m3sig, par_p3sig))
    return


def calcBestModel(samples, freeParID, fixedParams,
                  wlObs, specIobs, errObs):
    """
    Run and save a model with the final 'best' parameter values.
    Calculate the best parameters from the median of the good (after burn-in)
    portion of the chain (in the samples array).
    Then run another final model with lnlike, using verbose to save outout.
    """
    finalParams = np.percentile(samples, 50., axis=0)
    lnlikeEnd = lnlike(finalParams, freeParID, fixedParams,
                       wlObs, specIobs, errObs, verbose=True)
    return lnlikeEnd


def runMCMC(nchain, nburn, nwalkers, params, epsilons, fixedParams, wlRanges,
            obsName='observed.dat', verbose=False):
    """
    Run the main MCMC routine, calling emcee.

    This initiates and runs the MCMC chain, saves the results to chain.dat,
    prints median parameter values from the chain after burn-in,
    and saves a synthetic spectrum calculated with those median values
    (interpolated onto the observed wavelength pixels) to outSpeci.dat

    Input parameters use Python dictionaries with keyword names for model
    parameters.   Know dictionary keywords are:
    Vr (km/s);  vsini (km/s);  Vmic (km/s);  Vmac (km/s);
    Teff;  logg;  metal (metallicity);
    contFlx (additional veiling continuum flux);
    abun_## (## = atomic number);
    contN_# (# = polynomial degree, e.g. 0, 1, ...) for geometric polynomials;
    or alternatively contC_# for Chebyshev polynomials;
    or alternatively contL_# for Legendre polynomials;
    Bmono (provide a List of values); FFmono (list); Bdip (list); FFdip (list);
    errLogF the log of a scaling factor applied to the observed uncertainties
    used as: err = errObs*np.exp(errLogF)
    
    Takes:
    nchain - the total length of the MCMC chain, used by emcee
    nburn - the number of burn in steps in the MCMC chain, used by emcee
        (the number of steps discarded)
    nwalkers - the number of MCMC walkers used by emcee
    params - dictionary of free model parameters.  
        The dictionary has parameter key word names and values.
    epsilons - dictionary of 1sigma values for the initial distribution
        of parameters, used to initialize emcee's walkers.
        There must be an epsilon value for each free parameter in params.
        The dictionary has parameter key word names and values.
    fixedParams - dictionary of model parameters held fixed.
        The dictionary has parameter key word names and values.
    wlRanges - list of lists, containing sets of start and end wavelength
        values for the portions of spectrum fit.  Multiple wavelength ranges
        can be supplied.  Wavelengths are in the same rest frame as
         the observation.
    obsName - name of the observed spectrum file
    """

    freePar, freeParID, freeParEps = parseInputParams(params, epsilons, fixedParams)

    # Read in the observation, and get the pixels in the range(s) to be fit
    wlObs, specIobs, errObs = readObs.getObsSpecInRange(wlRanges, obsName)
    
    # Setup input for emcee
    # Set random initial walker positions for emcee
    npars = len(freePar)
    pos0 = [freePar + freeParEps*np.random.randn(npars) for i in range(nwalkers)]

    # Initialze the MCMC sampler
    sampler = emcee.EnsembleSampler(nwalkers, npars, lnprob,
                    args=(freeParID, fixedParams, wlObs, specIobs, errObs))
    # Setup the output file
    writeChainFileHeader(freeParID, fixedParams, wlRanges)
    
    #Full MCMC run, iterate through saving results
    iStep = 0
    for state in sampler.sample(pos0, iterations=nchain, store=True):
        print('step {:}'.format(iStep+1))
        position, lnprobability, rng = state
        writeChainIteration(position, lnprobability)
        
        iStep += 1

    # Cut burn-in part of the chain, and flatten it over walkers
    if nburn > nchain or nburn < 0:
        print('ERROR: bad nburn value! {:}'.format(nburn))
        nburn = 0
    print('burning {:} steps ({:} samples)'.format(nburn, nburn*nwalkers))
    samples = sampler.chain[:, nburn:, :].reshape((-1, npars))

    # Generate a final 'best' model
    calcBestModel(samples, freeParID, fixedParams, wlObs, specIobs, errObs)
    
    # Print out ending information and the final best parameters
    printFinalParams(samples, freeParID, fixedParams,
                     lnprobability, sampler=sampler, verbose=False)
    
    return sampler.chain


def continueMCMC(chainFile, nadd, nburn,
                 obsName='observed.dat', verbose=False):
    """
    Add additional steps to an existing MCMC chain.
    This reads a chain file, then calculates more steps by calling emcee.
    """
    
    # Read the chain file, and get both the current parameter values
    # and parameter names.
    freeParID, samples0, npars, nsteps0, nwalkers, fixedParams, wlRanges = readChain(chainFile)
    
    # Read in the observation, and get the pixels in the range(s) to be fit
    wlObs, specIobs, errObs = readObs.getObsSpecInRange(wlRanges, obsName)
    
    # Setup input for emcee
    # Set initial walker positions from end of current chain
    pos0 = samples0[-nwalkers:, :]
    
    # Initialze the MCMC sampler
    sampler = emcee.EnsembleSampler(nwalkers, npars, lnprob,
                    args=(freeParID, fixedParams, wlObs, specIobs, errObs))
    
    #Full MCMC run, iterate through saving results
    iStep = 0
    for state in sampler.sample(pos0, iterations=nadd, store=True):
        print('step {:}'.format(iStep+1))
        position, lnprobability, rng = state
        writeChainIteration(position, lnprobability, fileName=chainFile)
        
        iStep += 1

    # Cut burn-in part of the chain, and flatten it over walkers
    if nburn > nadd + nsteps0 or nburn < 0:
        print('ERROR: bad nburn value! {:}'.format(nburn))
        nburn = 0
    print('burning {:} steps ({:} samples)'.format(nburn, nburn*nwalkers))
    samples = sampler.chain.reshape((-1, npars))
    samples = np.concatenate((samples0, samples), axis=0)[nburn*nwalkers:, :]
    
    # Generate a final 'best' model
    calcBestModel(samples, freeParID, fixedParams, wlObs, specIobs, errObs)
    
    # Print out ending information and the final best parameters
    printFinalParams(samples, freeParID, fixedParams,
                     lnprobability, sampler=sampler, verbose=False)

    return sampler.chain
