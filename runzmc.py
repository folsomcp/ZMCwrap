import zmcwrap

#### stellar parameters ####
# Know dictionary keywords are:
# Vr,  vsini,  Vmic,  Vmac,  Teff,  logg,  metal,  contFlx,
# abun_## (## = atomic number),
# contN_# (# = polynomial degree, e.g. 0, 1, ...) for geometric polynomials,
# or contC_# for Chebyshev polynomials,
# or contL_# for Legendre polynomials,
# Bmono (provide a List), FFmono (List), Bdip (List), FFdip (List),
# errLogF the log of a scaling factor applied to the observed uncertainties
#     used as: err = errObs*np.exp(errF)

# params defines initial parameters passed to the fitting routine.  
# Omitted values will be read from the input file (inlmam.dat).
# epsilons defines the initial (Gaussian) distribution of parameters
# for the ensemble of walkers in the MCMC routine.
# For a parameter to be 'fitable' (free in this modelling) it needs a
# value in the params dict and in the epsilons dict.
# For a parameter to be 'fixed' it can be set it in the fixedParams dict,
# or it can be omitted here and set in the usual Zeeman input files.

params =   {'Teff':9000.0, 'Vr':5.0, 'vsini':30.0, 'Vmic':2.0, 
            'abun_2': -1.07,
            'abun_6': -3.57,
            'abun_8': -3.31,
            'abun_11':-5.76,
            'abun_12':-4.40,
            'abun_14':-4.49,
            'abun_16':-4.88,
            'abun_20':-5.66,
            'abun_21':-8.85,
            'abun_22':-7.05,
            'abun_24':-6.36,
            'abun_26':-4.50,
            'abun_28':-5.78,
            'errLogF':0.0,
            'contC_0':1.0, 'contC_1':0.0, 'contC_2':0.0, 'contC_3':0.0}
epsilons = {'Teff':10., 'Vr':0.1, 'vsini':0.1, 'Vmic':0.1, 
            'abun_2': 0.01,
            'abun_6': 0.01,
            'abun_8': 0.01,
            'abun_11':0.01,
            'abun_12':0.01,
            'abun_14':0.01,
            'abun_16':0.01,
            'abun_20':0.01,
            'abun_21':0.01,
            'abun_22':0.01,
            'abun_24':0.01,
            'abun_26':0.01,
            'abun_28':0.01,
            'errLogF':0.01,
            'contC_0':1e-4, 'contC_1':1e-4, 'contC_2':1e-4, 'contC_3':1e-4}
fixedParams = {'logg':4.3, 'Vmac':0.0,
               'abun_56':-9.82,}

wlRanges = [[5390.0, 5880.0], [5955.0, 6000.0]]

nchain = 500
nburn = 200
nwalkers = 200

chain = zmcwrap.runMCMC(nchain, nburn, nwalkers, 
                        params, epsilons, fixedParams, wlRanges,
                        obsName='observed.dat')

## Optionally, generate model spectra for the last n samples in the chain
#saveLastSpec = 1000
#zmcwrap.calcSynSpecSamples(saveLastSpec, nthin=10,
#                           chainFile='chain.dat', obsFile='observed.dat',
#                           specSamplesFile='outSpecSamples.dat')
