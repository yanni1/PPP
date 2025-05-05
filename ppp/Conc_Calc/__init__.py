# -*- coding: utf-8 -*-

"""
## Python (sub)module Conc_Calc
"""

import numpy as np
import params as p
from reac_rate import reaction_rate

#conc grad calc 
def update_conc(C_CO, C_CO2, C_O2, O2_reservoir, v):
    
    Cn_CO = np.copy(C_CO) #new array for next time step
    Cn_CO2 =np.copy(C_CO2)
    Cn_O2 = np.copy(C_O2)
    
    for j in range(1,p.ny-1): #exclude exterior points for now
        for i in range(1,p.nx-1):
            
            #diffusion
            D2x_CO = (C_CO[j, i+1] + C_CO[j, i-1] - 2 * C_CO[j, i]) / p.dx**2
            D2y_CO = (C_CO[j+1, i] + C_CO[j-1, i] - 2 * C_CO[j, i]) / p.dy**2
            Diff_CO = p.D*(D2x_CO + D2y_CO)
            
            D2x_CO2 = (C_CO2[j, i+1] + C_CO2[j, i-1] - 2 * C_CO2[j, i]) / p.dx**2
            D2y_CO2 = (C_CO2[j+1, i] + C_CO2[j-1, i] - 2 * C_CO2[j, i]) / p.dy**2
            Diff_CO2 = p.D*(D2x_CO2 + D2y_CO2)
            
            D2x_O2 = (C_O2[j, i+1] + C_O2[j, i-1] - 2 * C_O2[j, i]) / p.dx**2
            D2y_O2 = (C_O2[j+1, i] + C_O2[j-1, i] - 2 * C_O2[j, i]) / p.dy**2
            Diff_O2 = p.D*(D2x_O2 + D2y_O2)
                   
            #advection
            dy_CO = (C_CO[j, i] - C_CO[j-1, i]) / p.dy
            advec_CO = v[i] * dy_CO
            
            dy_CO2 = (C_CO2[j, i] - C_CO2[j-1, i]) / p.dy
            advec_CO2 = v[i] * dy_CO2
            
            dy_O2 = (C_O2[j, i] - C_O2[j-1, i]) / p.dy
            advec_O2 = v[i] * dy_O2
            
            #reaction
            
            k = reaction_rate(C_CO[j,i], C_CO2[j,i], C_O2[j,i])
            reac_CO = k * C_CO[j,i]
            reac_CO2 = -k * C_CO[j,i] #right side of reaction
            reac_O2 = k * C_CO[j,i]
            
            #total
            Cn_CO[j,i] = C_CO[j,i] + (p.dt * (Diff_CO - advec_CO - reac_CO))
            Cn_CO2[j,i] = C_CO2[j,i] + (p.dt * (Diff_CO2 - advec_CO2 - reac_CO2))
            Cn_O2[j,i] = C_O2[j,i] + (p.dt * (Diff_O2 - advec_O2 - reac_O2))#should this be times p.dt or not?
            #O2 circle
            dist_sq = ((i - p.a) ** 2) + ((j - p.b) ** 2) #distance from the center of the circle
            if dist_sq <= p.r2: #check if point is inside circle
                Cn_O2[j,i] = O2_reservoir
    #edges
    Cn_CO[:,0] = Cn_CO[:, 1]
    Cn_CO[:, -1] = Cn_CO[:, -2]
    Cn_O2[:, 0] = Cn_O2[:, 1]
    Cn_O2[:, -1] = Cn_O2[:, -2]
    Cn_CO2[:, 0] = Cn_CO2[:, 1]
    Cn_CO2[:, -1] = Cn_CO2[:, -2]
    
   
    return Cn_CO, Cn_CO2, Cn_O2
