#!/usr/bin/env python3
import os
from pathlib import Path
import sys

__projectdir__ = Path(os.path.dirname(os.path.realpath(__file__)) + '/')

import matplotlib.pyplot as plt
import numpy as np

# Vary Parameters Under Two Sector Model:{{{1
def varyrigidityonesector():
    rigiditylist = np.linspace(0.05, 0.5, 100)
    aggMClist = []
    for rigidity2 in rigiditylist:
        from manysector_ss_func import getss
        retdict = getss(0.96 **(1/12), [0.1, rigidity2], 1.02**(1/12), 8, 8, [0.5, 0.5])
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
        from manysector_ss_func import getss
        retdict = getss(0.96 **(1/12), [0.2, 0.1], 1.02**(1/12), 8, tau, [0.4, 0.6])
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
    sys.path.append(str(__projectdir__ / Path('submodules/calvo-basic-ss/')))
    from calvo_ss_func import calvobasicss
    _, MC, NU = calvobasicss(BETA, LAMBDA, SIGMA, 1)
    fbasic = 1 - MC * NU - profitsharestar
    finter = 1 - MC * NU - profitsharestar * (1 - s_m * MC * NU)

    profitsbasic = []
    profitsinter = []
    for Pistar in Pistars_monthly:
        from manysector_ss_func import calvobasicss
        _, MC, NU = calvobasicss(BETA, LAMBDA, SIGMA, Pistar)
        profitsbasic.append(1 - MC * NU - fbasic)
        profitsinter.append( (1 - MC * NU - finter) / (1 - s_m * MC * NU) )

    # get fixed costs for multisector calvo 
    from manysector_ss_func import ns_ss
    retdict = ns_ss(BETA, 1, SIGMA, SIGMA)
    aggMC = retdict['MC'] * retdict['NU']
    fmulti = 1 - aggMC - profitsharestar
    fmultiinter = 1 - aggMC - profitsharestar * (1 - s_m * MC * NU)

    profitsmulti = []
    profitsmultiinter = []
    for Pistar in Pistars_monthly:
        from manysector_ss_func import ns_ss
        retdict = ns_ss(BETA, Pistar, SIGMA, SIGMA)
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
