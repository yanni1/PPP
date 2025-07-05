# Animation update function
import Calc_CPP # pyright: ignore[reportMissingImports]
import params_module # pyright: ignore[reportMissingImports]
import Genetic_algo # pyright: ignore[reportMissingImports]
import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation
from matplotlib.patches import Ellipse
import numpy as np
from calc_starter import get_best_params
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


nt = 5000

tau, rho = get_best_params(nt)
if rank == 0:
    print("GA found best tau, rho to be: ",tau,rho)
    p = params_module.params(nt, int(tau), rho)

    C_flat, eps_field = Calc_CPP.Calc_CPP(p)

    # C needs to be a 4D array: C[t][s][y][x]
    C = np.array(C_flat, dtype=np.float32).reshape((p.nt, p.ns, p.ny, p.nx))
    eps_field = np.array(eps_field  , dtype=np.float32).reshape((p.nt, p.ny, p.nx))

    extent = [0, p.nx, 0, p.ny]

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4)
    ax1.set_title("CO")
    ax2.set_title("CO₂")
    ax3.set_title("O₂")
    ax4.set_title("ε_CO (absorption)")
    ims = []

    # ellips
    x_center = p.x0_cc
    y_center = p.y0_cc
    width = 2 * p.semiMaj
    height = 2 * p.semiMin

    #ellipse patches
    ell1 = Ellipse((x_center, y_center), width, height, edgecolor='orange', facecolor='none', lw=1.5, linestyle='--')
    ell2 = Ellipse((x_center, y_center), width, height, edgecolor='orange', facecolor='none', lw=1.5, linestyle='--')
    ax1.add_patch(ell1)
    ax2.add_patch(ell2)

    #loop over time steps
    for t in range(nt):
        if t % 100 == 0:
            con = C[t]
            eps = eps_field[t]
            # con is a 3D array at this time step: con[s][y][x]
            im1 = ax1.imshow(con[0], animated=True, origin='lower', extent=extent, aspect='equal')
            im2 = ax2.imshow(con[1], animated=True, origin='lower', extent=extent, aspect='equal')
            im3 = ax3.imshow(con[2], animated=True, origin='lower', extent=extent, aspect='equal')
            im4 = ax4.imshow(eps, cmap='inferno', origin='lower', extent=extent, aspect='equal')
            ims.append([im1, im2, im3, im4])
        else: 
            continue

    fig.colorbar(ax4.images[0], ax=ax4, orientation='vertical', label='Concentration')

    ani = ArtistAnimation(fig, ims, 100)
    plt.show()


# mpirun -np 2 /Users/yanni/.pyenv/versions/3.12.2/bin/python /Users/yanni/VSC_Projects_Folder/testenv/PPP/ppp/Animation_CPP.py