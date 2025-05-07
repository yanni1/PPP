# -*- coding: utf-8 -*-

"""
## Python (sub)module params
"""

import numpy as np

# Parameters
dx = 0.01 #m
dy = 0.01 #m
dt = 0.0005 #s
Lx = 0.5 #m
Ly = 2.0 #m
D = 0.002 #m^2/s
k0 = 1000 #s^-1
vmax = 0.5 #m/s
CO_reservoir = 1.2 #mol/m^3
O2_reservoir = 1.0 #mol/m^3

#Grid
nx = int(Lx/dx) #number of gridpoints in x
ny = int(Ly/dy) #number of gridpoints in y
nt =100 #5s time interval
x = np.linspace(0, Lx, nx) #x_grid array
y = np.linspace(0, Ly, ny) #y_grid array
#o2 vent coords
b = ny/10 #center of cicle y
a = nx/4 #center of circle x
r = nx/8 #radius circle
r2 = r ** 2 #radius squared
