#!/usr/bin/env python3
import os
from pathlib import Path
import sys

__projectdir__ = Path(os.path.dirname(os.path.realpath(__file__)) + '/')

import numpy as np

# Get Steady State:{{{1
def getss(BETA, LAMBDAs, Pistar, SIGMA, TAU, WEIGHTs):
    """
    Compute the steady state for a multiple sector model.

    SIGMA = elasticity of substitution within sectors, TAU = across sectors
    """

    J = len(WEIGHTs)

    # verify WEIGHTs sum to 1
    if np.abs(np.sum(WEIGHTs) - 1) > 1e-6:
        raise ValueError('WEIGHTs must sum to 1.')

    # compute PjstaroverPj
    PjstaroverPj_list = []
    for j in range(J):
        PstaroverPj = ( (1-(1-LAMBDAs[j])*Pistar**(SIGMA-1)) / LAMBDAs[j])**(1/(1-SIGMA))
        if not PstaroverPj > 0:
            raise ValueError(f"PstaroverPj = {PstaroverPj}. The model will not solve when this value is less than zero. This happens for sector j: {j}. It is because a small fraction of that sectors prices change each period ({LAMBDAs[j]}) while the remainder falls by Pistar ({Pistar}) in steady state. Pistar is large enough that PstaroverPj cannot be large enough to make these equalize. To prevent this, I need to lower Pistar, lower SIGMA, or raise LAMBDAj.")
        PjstaroverPj_list.append( PstaroverPj )

    # compute nu_j
    NUj_list = []
    for j in range(J):
        NUj_list.append( 1/(1 - (1-LAMBDAs[j])*Pistar**SIGMA) * LAMBDAs[j] * PjstaroverPj_list[j]**(-SIGMA) )

    # compute MC
    # first compute the sum of the integral term to (1 - TAU)
    sumterm = 0
    terminfoc_list = []
    for j in range(J):
        terminfoc_part1 = SIGMA / (SIGMA - 1) * (1 - (1-LAMBDAs[j])*BETA*Pistar**(SIGMA-1)) / ( 1 - (1-LAMBDAs[j]) * BETA * Pistar**SIGMA )
        terminfoc_part2 = 1 / PjstaroverPj_list[j] 
        terminfoc = terminfoc_part1 * terminfoc_part2
        terminfoc_list.append(terminfoc)
        if terminfoc < 0:
            raise ValueError('Solving for MC failed in Calvo multisector. terminfoc < 0.')

        sumterm = sumterm + WEIGHTs[j] * (terminfoc_list[j]) ** (1 - TAU)
    if sumterm <= 0:
        raise ValueError('Solving for MC failed in Calvo multisector. sumterm < 0')
    MC = (1/sumterm)**(1/(1-TAU))

    # PjoverP_list
    PjoverP_list = []
    for j in range(J):
        PjoverP_list.append( terminfoc_list[j] * MC )

    # NU
    NU = 0
    for j in range(J):
        NU = NU + WEIGHTs[j] * PjoverP_list[j] ** (-TAU) * NUj_list[j]

    retdict = {}
    retdict['NUj_list'] = NUj_list
    retdict['PjstaroverPj_list'] = PjstaroverPj_list
    retdict['MC'] = MC
    retdict['PjoverP_list'] = PjoverP_list
    retdict['NU'] = NU

    if np.any([NU_j < 0 for NU_j in retdict['NUj_list']]):
        raise ValueError('NU_j take negative values.')

    return(retdict)


def test():
    LAMBDAs = [1 - (1 - 0.6) ** (1/4)]
    # LAMBDAs = [1 - (1 - 0.6) ** (1/4), 0.9]

    WEIGHTs = [1]

    retdict = getss(BETA = 0.94 ** (1/4), LAMBDAs = LAMBDAs, Pistar = 1.04 ** (1/4), SIGMA = 8, TAU = 1.001, WEIGHTs = WEIGHTs)
    print(retdict['MC'])
    print(retdict['PjoverP_list'])
    print(retdict['NUj_list'])
    print(retdict['NU'] * retdict['MC'])


def test2():
    LAMBDAs = [1 - (1 - 0.6) ** (1/4)] * 2
    # LAMBDAs = [1 - (1 - 0.6) ** (1/4), 0.9]

    WEIGHTs = [0.5, 0.5]

    retdict = getss(BETA = 0.94 ** (1/4), LAMBDAs = LAMBDAs, Pistar = 1.04 ** (1/4), SIGMA = 8, TAU = 1.001, WEIGHTs = WEIGHTs)
    print(retdict['MC'])
    print(retdict['PjoverP_list'])
    print(retdict['NUj_list'])
    print(retdict['NU'] * retdict['MC'])


# Pricing Parameters Based upon Nakamura Steinsson 2008:{{{1
def ns_vectors(numsectors = 14):
    """
    Vectors for weights and lambdas from 
    """
    if numsectors == 6:
        # Table 2, p.978 2010 multisector menu cost
        ns_weights = np.array([7.7, 19.1, 5.9, 13.7, 38.5, 15.1])
        ns_monthlyfreqs = np.array([91.6, 35.5, 25.4, 11.9, 8.8, 5.2])
    elif numsectors == 9:
        # Table 2, p.978 2010 multisector menu cost
        ns_weights = np.array([7.7, 19.1, 5.9, 9.2, 13.7, 9.6, 10.0, 15.1, 9.7])
        ns_monthlyfreqs = np.array([91.6, 35.5, 25.4, 19.7, 11.9, 7.6, 5.5, 5.2, 3.2])
    elif numsectors == 11:
        # weights from Table 2, p.1433 of Nakamura Steinsson's Five Facts About Prices, 2008
        ns_weights = np.array([8.2, 5.9, 5.0, 6.5, 8.3, 3.6, 5.4, 5.3, 5.1, 5.5, 38.5])
        ns_monthlyfreqs = np.array([10.5, 25.0, 6.0, 3.6, 31.3, 6.0, 15.0, 38.1, 87.6, 41.7, 6.1])
    elif numsectors == 14:
        # Table 2, p.978 2010 multisector menu cost
        ns_weights = np.array([7.7, 5.3, 5.5, 5.9, 8.3, 7.7, 13.7, 7.5, 5.0, 7.8, 3.6, 7.6, 6.5, 7.9])
        ns_monthlyfreqs = np.array([91.6, 49.4, 43.7, 25.4, 21.3, 21.7, 11.9, 8.4, 6.5, 6.2, 6.1, 4.9, 3.6, 2.9])
    else:
        raise ValueError('Incorrect option for sectors.')

    # adjustments
    ns_weights = ns_weights / np.sum(ns_weights)
    ns_monthlyfreqs = ns_monthlyfreqs / 100

    return(ns_weights, ns_monthlyfreqs)


def getns_lambdas(monthsinperiod = 3, numsectors = 14):
    ns_weights, ns_monthlyfreqs = ns_vectors(numsectors = numsectors)
    ns_lambdas = 1 - (1 - ns_monthlyfreqs) ** monthsinperiod
    return(ns_weights, ns_lambdas)


def ns_ss(BETA, Pistar, SIGMA, TAU, monthsinperiod = 1):
    """
    Frequencies of price changes and weights taken from Nakamura Steinsson
    """
    ns_weights, ns_lambdas = getns_lambdas(monthsinperiod = monthsinperiod)

    retdict = getss(BETA, ns_lambdas, Pistar, SIGMA, TAU, ns_weights)

    return(retdict)
