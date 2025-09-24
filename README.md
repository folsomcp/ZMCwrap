# ZMCwrap

ZMCwrap is a tool for estimating stellar parameters using MCMC and spectrum synthesis.  It is built as a Python wrapper around the Zeeman spectrum synthesis code.

This requires the emcee package to run the MCMC sampler (see [emcee's docs](https://emcee.readthedocs.io/en/stable/)).  This also requires a copy of the Zeeman spectrum synthesis code to run, but could be adapted to other codes.

If you use this software in your research, please cite Folsom et al. (in prep.)

## Usage

There are a few useful functions provided by this package.  The main code is in the zmcwrap folder, while several executable scripts are in the root folder.

The `runzmc.py` file provides a example of setting up and executing an MCMC modelling run.  In this example the observation file is read from `observed.dat`.  The observation file should be a text file with columns of wavelength, flux, and uncertainties on the flux; or columns of wavelength, flux, three columns that are skipped (e.g. polarimetric information) and uncertainties on the flux.

Running this function requires having a version of Zeeman set up, because functions used by `runzmc.py` will attempt to run Zeeman.  Some important information must be set up with Zeeman.  Particularly the line list and any desired oscillator strength corrections must be set up as usual.  The wavelength range for the synthetic spectrum must also be set in Zeeman, and this should be somewhat wider than the wavelength range set for fitting in `runzmc.py` (with the `wlRanges` list).  The synthetic spectrum wavelength range should be somewhat wider to allow for Doppler shifting by a radial velocity.  

The `continuezmc.py` file provides an example of continuing an existing MCMC chain, and adding more samples to the end of the chain.  In this example, the script uses the existing `chain.dat` file.

The `plotChain.py` file plots parameter values as a function of step in the chain.  Then it prints out the median values, as well as uncertainties from percentiles corresponding to a 1sigma and 3sigma confidence level.  Then it generates a corner plot of the parameter distributions.  When running this script it is important to specify the number of burn-in steps in the chain.  This is done by running the script with the `-b` flag, like: `python plotChain.py -b 300`.  You can run the script with the `-h` flag (`python plotChain.py -h`) for more help information.

The `plotChainAltLabels.py` behaves like `plotChain.py` although it uses alternative shorter labels for the parameters.

The `plotSpecSamp.py` can be used if you have run the `zmcwrap.calcSynSpecSamples` function.  This will plot the observation, best fit (median) model, and the set of individual sample models.  You can run the script with `python plotSpecSamp.py -h` for some help information.  By default the script will try to read an observation from `observed.dat` or it can be specified with `-o`.  The default best fit model is `outSpeci.dat` or it can be specified with `- b`.  The default set of sample spectra is `outSpecSamples.dat` or it can be specified with `-s`.

The `plotSpecSampLines.py` is similar to `plotSpecSamp.py`, but it also can plot a line list.  This uses the SpecpolFlow package (see [SpecpolFlow's docs](https://folsomcp.github.io/specpolFlow/index.html)) to plot the line list.  The line list can be specified with `-l`, and use `-h` for more information.

## Main Functions

The ZMCwrap package provide a few main functions.  These are used in the executable scripts mentioned above.  `runMCMC` runs the main MCMC code. `continueMCMC` adds more samples to an existing MCMC chain.  `readChain` reads a chain file for further processing.  `chainPlot` and `cornerPlotLight` both plot parameters as a function of step in the chain, and then generate corner plots of the parameter distributions.  `calcSynSpecSamples` calculates a set of synthetic spectra corresponding to the last n steps in the chain, and `readSpecSamples` reads back in the sample spectra (returning them as numpy arrays).  

## Adapting

While ZMCwrap is designed to run with the Zeeman spectrum synthesis code, it can be modified to use other spectrum synthesis codes to generate model spectra.  Doing this requires replacing the current `runSpecSynth` function in the `zmcwrap/wrapper_zeeman.py` file.  This function needs to take lists of free and fixed parameters, and return lists with the wavelength and continuum normalized flux for the model spectrum (other functions will interpolate it onto the observed spectrum).  Other functions in that file can be used a an example of how to write input files used by the spectrum synthesis code.  