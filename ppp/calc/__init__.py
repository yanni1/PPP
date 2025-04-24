# -*- coding: utf-8 -*-

"""
## Python (sub)module calc
"""

import numpy as np
import time
from Conc_Calc import update_conc
import params as p

t0 = time.time()

#Velocity Porfile   
v = p.vmax * (1 - 0.5*(((p.x-(p.Lx/2))/p.Lx)**2))

nt = p.nt
def calc(nt=nt):
    t2 = time.time()
    C_O2 = np.zeros((p.ny,p.nx)) #make conc matrix
    C_CO = np.zeros((p.ny,p.nx))
    C_CO2 = np.zeros((p.ny,p.nx))
    
    #boudary cond Co at bottom
    C_CO[0,:] = p.CO_reservoir
    #O2 vent initial condition
    #center of circle 
    b = p.ny/10 #center of cicle y
    a = p.nx/4 #center of circle p.x
    r = p.nx/8 #radius circle
    r2 = r ** 2 #radius squared
    O2_source = np.zeros((p.ny,p.nx))

    C = []  
    for N in range(p.nt):
        C_CO, C_CO2, C_O2 = update_conc(C_CO, C_CO2, C_O2, p.O2_reservoir,v)
        C.append((C_CO,C_CO2,C_O2))
    t3 = time.time()
    tijdCalc = t3-t2
    print('runtime',tijdCalc)
    return C






