# -*- coding: utf-8 -*-

"""
## Python (sub)module params
"""

from dataclasses import dataclass
import numpy as np

@dataclass
class p:
    dx: float = 0.01
    dy: float = 0.01
    dt: float = 0.0005
    Lx: float = 0.5
    Ly: float = 2.0
    D: float = 0.002
    k0: float = 1000
    vmax: float = 0.5
    CO_reservoir: float = 1.2
    O2_reservoir: float = 1.0

    def __post_init__(self): #function members
        self.nx = int(self.Lx / self.dx)
        self.ny = int(self.Ly / self.dy)
        self.nt = 100
        self.x = np.linspace(0, self.Lx, self.nx)
        self.y = np.linspace(0, self.Ly, self.ny)
        self.a = self.nx * 0.25
        self.b = self.ny * 0.1
        self.r = self.nx * 0.125
        self.r2 = self.r * self.r
