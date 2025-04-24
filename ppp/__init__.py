# -*- coding: utf-8 -*-

"""
## Python package ppp
"""

__version__ = '0.0.0'

from params import k0
#calc k
def reaction_rate(C_CO, C_CO2, C_O2):
    k_b = k0/4 * C_CO2
    k_f = k0 * C_CO * (C_O2 ** 2)
    k_consumption = k_f - k_b
    print(type(k_consumption))
    return k_consumption
