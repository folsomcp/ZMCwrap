import zmcwrap


# Number of steps to add to the existing chain chain,
# and steps to burn from the full final chain
nadd = 500
nburn = 500

chain = zmcwrap.continueMCMC('chain.dat', nadd, nburn)

## Optionally, generate model spectra for the last n samples in the chain
#saveLastSpec = 1000
#zmcwrap.calcSynSpecSamples(saveLastSpec, nthin=10,
#                           chainFile='chain.dat', obsFile='observed.dat',
#                           specSamplesFile='outSpecSamples.dat')
