# -*- coding: utf-8 -*-

"""
## Python (sub)module reac_rate
"""

from params import k0
#calc k
def reaction_rate(C_CO, C_CO2, C_O2):
    k_b = k0/4 * C_CO2
    k_f = k0 * C_CO * (C_O2 ** 2)
    k_consumption = k_f - k_b
    return k_consumption

