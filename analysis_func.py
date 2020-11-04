#!/usr/bin/env python3
# PYTHON_PREAMBLE_START_STANDARD:{{{

# Christopher David Cotton (c)
# http://www.cdcotton.com

# modules needed for preamble
import importlib
import os
from pathlib import Path
import sys

# Get full real filename
__fullrealfile__ = os.path.abspath(__file__)

# Function to get git directory containing this file
def getprojectdir(filename):
    curlevel = filename
    while curlevel is not '/':
        curlevel = os.path.dirname(curlevel)
        if os.path.exists(curlevel + '/.git/'):
            return(curlevel + '/')
    return(None)

# Directory of project
__projectdir__ = Path(getprojectdir(__fullrealfile__))

# Function to call functions from files by their absolute path.
# Imports modules if they've not already been imported
# First argument is filename, second is function name, third is dictionary containing loaded modules.
modulesdict = {}
def importattr(modulefilename, func, modulesdict = modulesdict):
    # get modulefilename as string to prevent problems in <= python3.5 with pathlib -> os
    modulefilename = str(modulefilename)
    # if function in this file
    if modulefilename == __fullrealfile__:
        return(eval(func))
    else:
        # add file to moduledict if not there already
        if modulefilename not in modulesdict:
            # check filename exists
            if not os.path.isfile(modulefilename):
                raise Exception('Module not exists: ' + modulefilename + '. Function: ' + func + '. Filename called from: ' + __fullrealfile__ + '.')
            # add directory to path
            sys.path.append(os.path.dirname(modulefilename))
            # actually add module to moduledict
            modulesdict[modulefilename] = importlib.import_module(''.join(os.path.basename(modulefilename).split('.')[: -1]))

        # get the actual function from the file and return it
        return(getattr(modulesdict[modulefilename], func))

# PYTHON_PREAMBLE_END:}}}

import matplotlib.pyplot as plt
import numpy as np

# Vary Parameters Under Two Sector Model:{{{1
def varyrigidityonesector():
    rigiditylist = np.linspace(0.05, 0.5, 100)
    aggMClist = []
    for rigidity2 in rigiditylist:
        retdict = importattr(__projectdir__ / Path('manysector-ss_func.py'), 'getss')(0.96 **(1/12), [0.1, rigidity2], 1.02**(1/12), 8, 8, [0.5, 0.5])
        aggMClist.append(1 - retdict['MC'] * retdict['NU'])

    plt.plot(rigiditylist, aggMClist)
    plt.xlabel('Frequency of Price Change in Second Sector')
    plt.ylabel('Profit Share')

    plt.show()

    plt.savefig(__projectdir__ / Path('temp/profitshare_frequencyonesector.png'))

    plt.clf()


def varytau():
    """
    Changing tau doesn't seem to make any difference.
    """
    taulist = np.linspace(2, 10, 100)
    aggMClist = []
    for tau in taulist:
        retdict = importattr(__projectdir__ / Path('manysector-ss_func.py'), 'getss')(0.96 **(1/12), [0.2, 0.1], 1.02**(1/12), 8, tau, [0.4, 0.6])
        aggMClist.append(1 - retdict['MC'] * retdict['NU'])

    plt.plot(taulist, aggMClist)
    plt.xlabel('TAU')
    plt.ylabel('Profit Share')

    plt.show()

    plt.savefig(__projectdir__ / Path('temp/tau_profitshare.png'))

    plt.clf()


# Inflation and Profit Share Different Calvo Models:{{{1
def inflation_profitshare_inter():
    Pistars = np.linspace(0.99, 1.04, 21)
    Pistars_monthly = Pistars ** (1/12)

    # profit share at zero inflation that I target
    profitsharestar = 0.1
    BETA = 0.96 ** (1/12)
    SIGMA = 8
    
    LAMBDA = 0.087
    s_m = 8/7 * 0.52
    
    # get fixed costs for basic calvo 
    _, MC, NU = importattr(__projectdir__ / Path('submodules/calvo-basic-ss/calvo-ss_func.py'), 'calvobasicss')(BETA, LAMBDA, SIGMA, 1)
    fbasic = 1 - MC * NU - profitsharestar
    finter = 1 - MC * NU - profitsharestar * (1 - s_m * MC * NU)

    profitsbasic = []
    profitsinter = []
    for Pistar in Pistars_monthly:
        _, MC, NU = importattr(__projectdir__ / Path('manysector-ss_func.py'), 'calvobasicss')(BETA, LAMBDA, SIGMA, Pistar)
        profitsbasic.append(1 - MC * NU - fbasic)
        profitsinter.append( (1 - MC * NU - finter) / (1 - s_m * MC * NU) )

    # get fixed costs for multisector calvo 
    retdict = importattr(__projectdir__ / Path('manysector-ss_func.py'), 'ns_ss')(BETA, 1, SIGMA, SIGMA)
    aggMC = retdict['MC'] * retdict['NU']
    fmulti = 1 - aggMC - profitsharestar
    fmultiinter = 1 - aggMC - profitsharestar * (1 - s_m * MC * NU)

    profitsmulti = []
    profitsmultiinter = []
    for Pistar in Pistars_monthly:
        retdict = importattr(__projectdir__ / Path('manysector-ss_func.py'), 'ns_ss')(BETA, Pistar, SIGMA, SIGMA)
        aggMC = retdict['MC'] * retdict['NU']
        profitsmulti.append(1 - aggMC - fmulti)
        profitsmultiinter.append( (1 - aggMC - fmultiinter) / (1 - s_m * aggMC) )


    plt.plot(Pistars, profitsbasic, label = 'Basic')
    plt.plot(Pistars, profitsinter, label = 'Intermediate Goods')
    plt.plot(Pistars, profitsmulti, label = 'Multisectors')
    plt.plot(Pistars, profitsmultiinter, label = 'Multisectors, Intermediate Goods')

    plt.legend()

    plt.savefig(__projectdir__ / Path('temp/multi_inflation_profitshare.png'))

    plt.show()

    plt.clf()


# Full:{{{1
def full():
    varyrigidityonesector()
    varytau()


# Run:{{{1
# full()
