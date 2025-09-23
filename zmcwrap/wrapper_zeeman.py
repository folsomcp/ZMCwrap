"""
Interface functions that provide a wrapper around the Zeeman spectrum
synthesis code.  
This could be adapted for other spectrum synthesis codes,
or even grids of pre-computed spectra.
"""

import numpy as np
import subprocess

zeemanCommand = './lma-sub'
#zeemanCommand = './lmamp'
zeemanInput = 'inlmam.dat'

#Modify the main input file for Zeeman, inlmam.dat
#Read the file into memory, modify any parameters,
#then write a new version back to disk.
def updateInlma(freePar, freeParID, fixedParams):
    """
    Read in and then re-write the inlmam.dat file.
    This takes the free parameters (and their keywords) from freePar and 
    freeParID, and fixed parameters (as a dict) in fixedParams.
    The inlmam.dat file must already exist and be in the correct format
    for Zeeman. All parameters will be set to fixed.
    Parameter values from inlmam.dat will be reused unless they are specified
    by the input function arguments.
    """
    #Read the Zeeman-LMA input
    fInlma = open(zeemanInput, 'r')
    
    txtInlma = []
    monoB = []
    monoff = []
    numMono = 0
    dipB = []
    dipff=[]
    numDip = 0
    numEl = 0
    listEl = []
    abuns = []
    stopEl = 0
    i=1
    for line in fInlma:
        if i == 4:
            vr = float(line.split()[0])
        if i == 6:
            vsini = float(line.split()[0])
        if i == 8:
            vmic = float(line.split()[0])
        if i == 10:
            vmac = float(line.split()[0])
        if i == 12:
            Teff = float(line.split()[0])
        if i == 14:
            logg = float(line.split()[0])
        if i == 16:
            metal = float(line.split()[0])
        if i == 18:
            contFlux = float(line.split()[0])
        if i == 20:
            numMono = int(line.split()[0])
        if (i > 21) and (i < 22+numMono):
            monoB += [float(line.split()[0])]
            monoff += [float(line.split()[1])]
        if i == 23 + numMono:
            numDip = int(line.split()[0])
        if (i > 24+numMono) and (i < 25+numMono+numDip):
            dipB += [float(line.split()[0])]
            dipff += [float(line.split()[1])]
        # (always set the continuum polynomial fitting in Zeeman to 0
        if (i > 25+2+numMono+numDip):
            if (stopEl == 0):
                if (line.strip()[0] == '='):
                    stopEl = 1
                else:
                    numEl += 1
                    listEl += [int(line.split()[0])]
                    abuns += [float(line.split()[1])]
        
        txtInlma += [line]
        i += 1
    
    fInlma.close()
    if(stopEl == 0): print('ERROR missing end of inlmam.dat (missing ending =)')

    
    #Set the new/modified parameters, overwriting the values read from file.
    if(len(freePar) != len(freeParID)):
        raise ValueError ('ERROR: missmatch in freePar and freeParID lengths')
    for key in freeParID + list(fixedParams):
        if key[:5] == 'abun_':
            abuns = []
            listEl = []
        elif key[:5] == 'Bmono': monoB = []
        elif key[:6] == 'FFmono': monoff = []
        elif key[:4] == 'Bdip': dipB = []
        elif key[:5] == 'FFdip': dipff = []
    
    for i in range(len(freePar)):
        if freeParID[i] == 'Vr':      vr       = freePar[i]*1e5
        if freeParID[i] == 'vsini':   vsini    = freePar[i]*1e5
        if freeParID[i] == 'Vmic':    vmic     = freePar[i]*1e5
        if freeParID[i] == 'Vmac':    vmac     = freePar[i]*1e5
        if freeParID[i] == 'Teff':    Teff     = freePar[i]
        if freeParID[i] == 'logg':    logg     = freePar[i]
        if freeParID[i] == 'metal':   metal    = freePar[i]
        if freeParID[i] == 'contFlx': contFlux = freePar[i]
        if freeParID[i][:5] == 'Bmono':  monoB  += [freePar[i]]
        if freeParID[i][:6] == 'FFmono': monoff += [freePar[i]]
        if freeParID[i][:4] == 'Bdip':   dipB   += [freePar[i]]
        if freeParID[i][:5] == 'FFdip':  dipff  += [freePar[i]]
        if freeParID[i][:5] == 'abun_':
            abuns   += [freePar[i]]
            listEl  += [int(freeParID[i][5:])]
    
    for fixID, fixPar in fixedParams.items():
        if fixID == 'Vr':      vr       = fixPar*1e5
        if fixID == 'vsini':   vsini    = fixPar*1e5
        if fixID == 'Vmic':    vmic     = fixPar*1e5
        if fixID == 'Vmac':    vmac     = fixPar*1e5
        if fixID == 'Teff':    Teff     = fixPar
        if fixID == 'logg':    logg     = fixPar
        if fixID == 'metal':   metal    = fixPar
        if fixID == 'contFlx': contFlux = fixPar
        if fixID[:5] == 'Bmono':  monoB  = fixPar
        if fixID[:6] == 'FFmono': monoff = fixPar
        if fixID[:4] == 'Bdip':   dipB   = fixPar
        if fixID[:5] == 'FFdip':  dipff  = fixPar
        if fixID[:5] == 'abun_':
            abuns   += [fixPar]
            listEl  += [int(fixID[5:])]
    
    #Extra error trapping (make sure things are lists and have the right length)
    if isinstance(monoB, (int, float)):  monoB  = [monoB]
    if isinstance(monoff, (int, float)): monoff = [monoff]
    if isinstance(dipB, (int, float)):   dipB   = [dipB]
    if isinstance(dipff, (int, float)):  dipff  = [dipff]
    if isinstance(listEl, (int, float)): listEl = [listEl]
    if isinstance(abuns, (int, float)):  abuns  = [abuns]
    if(len(monoB) != len(monoff)):
        raise ValueError ('ERROR inconsistent Bmono and FFmono lists: '
                          '{:} {:}'.format(monoB, monoff))
    if(len(dipB) != len(dipff)):
        raise ValueError ('ERROR inconsistent Bdip and FFdip lists: '
                          '{:} {:}'.format(dipB, dipff))
    if(len(listEl) != len(abuns)):
        raise ValueError ('ERROR inconsistent element and abundance lists:'
                          '{:} {:}'.format(listEl, abuns))
    
    #Write the Zeeman-LMA input
    fInlma = open(zeemanInput, 'w')
    i = 1
    for line in txtInlma:
        
        if i == 4:
            fInlma.write('{:14.8e}  {:1n}\n'.format(vr, 0))
        elif i == 6:
            fInlma.write('{:14.8e}  {:1n}\n'.format(vsini, 0))
        elif i == 8:
            fInlma.write('{:14.8e}  {:1n}\n'.format(vmic, 0))
        elif i == 10:
            fInlma.write('{:14.8e}  {:1n}\n'.format(vmac, 0))
        elif i == 12:
            fInlma.write('{:10.4f}      {:1n}\n'.format(Teff, 0))
        elif i == 14:
            fInlma.write('{:7.5f}         {:1n}\n'.format(logg, 0))
        elif i == 16:
            fInlma.write('{:8.5f}        {:1n}\n'.format(metal, 0))
        elif i == 18:
            fInlma.write('{:8.5f}        {:1n}\n'.format(contFlux, 0))
        elif i == 20:
            fInlma.write('{:2n}\n'.format(len(monoB)))
        elif i == 21:
            #Write the info header line and the new Bmono lines
            fInlma.write(line)
            for k in range(len(monoB)):
                fInlma.write('{:9.3f}    {:16.10e}   {:1n}       {:1n}\n'.format(
                    monoB[k], monoff[k], 0, 0))
        elif (i > 21) and (i < 22+numMono):
            #Do nothing for the old Bmono list
            pass
        elif i == 23 + numMono:
            fInlma.write('{:2n}\n'.format(len(dipB)))
        elif i == 24 + numMono:
            #Write the info header line and the new Bdip lines
            fInlma.write(line)
            for k in range(len(dipB)):
                fInlma.write('{:9.3f}    {:16.10e}   {:1n}       {:1n}\n'.format(
                    dipB[k], dipff[k], 0, 0))
        elif (i > 24+numMono) and (i < 25+numMono+numDip):
            #Do nothing for the old Bdip list
            pass
        elif i == 26+numMono+numDip:
            fInlma.write('{:<2n}\n'.format(0))
        elif i == 25+2+numMono+numDip:
            #Write the info header line and the new abundance lines
            fInlma.write(line)
            for k in range(len(listEl)):
                fInlma.write('{:<2n}       {:8.5f}  {:1n}\n'.format(
                    listEl[k], abuns[k], 0))
        elif (i > 25+2+numMono+numDip) and (i < 26+2+numMono+numDip+numEl):
            pass
        else:
            fInlma.write(line)
            
        i += 1
    fInlma.close()
    return


def getZeemanSpec(modelName = 'plotff1'):
    """
    Read in the output spectrum from Zeeman.
    Apply a Doppler shift here.  Interpolating onto a grid of observed points
    should be in the calling routine.
    """
    # With numpy version 2.0, try benchmarking
    #wl, specI = np.loadtxt(modelName, usecols=(0,1), unpack=True)
    wl = []
    specI = []
    fIn = open(modelName, 'r')
    for line in fIn:
        vals = line.split()
        wl.append(float(vals[0]))
        specI.append(float(vals[1]))
    fIn.close()
    
    wla = np.array(wl)
    # Before Zeeman version 9.7.17 this was needed (now it is included in plotff1)
    ## Apply the Doppler shift correction to the synthetic spectrum (in km/s)
    ##wla = wla + wla*vr/2.99792458e5
    return wla, np.array(specI)


def runSpecSynth(freePar, freeParID, fixedParams, verbose=False):
    """
    Run the spectrum synthesis code and return the resulting synthetic spectrum.
    This is built for Zeeman, specifically the lma version of the code.

    Takes:
    freePar - list of free parameter values
    freeParID - list of names for those free parameters
    fixedParams - dictionary of fixed parameter names and values
    verbose - a flag for printing extra diagnostic information

    Returns:
    wlSyn, specIsyn - arrays of wavelength and continuum normalized flux
    for the synthetic spectrum
    """

    #Calculate the model spectrum for this set of parameters with Zeeman
    updateInlma(freePar, freeParID, fixedParams)
    if verbose: print("Running Zeeman")
    subprocess.call(zeemanCommand)

    # for older versions of Zeeman that didn't Doppler shift plotff1
    ##Include a Doppler shift for radial velocity
    #if 'Vr' in freeParID:
    #    vr = freePar[freeParID.index('Vr')]
    #elif 'Vr' in fixedParams:
    #    vr = fixedParams['Vr']
    #else:
    #    vr = 0.0

    #Read the output from Zeeman
    wlSyn, specIsyn = getZeemanSpec()

    return wlSyn, specIsyn
