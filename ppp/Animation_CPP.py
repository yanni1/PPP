# Animation update function
import Calc_CPP
import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation

# C is a 4D array: C[t][s][y][x]
C = Calc_CPP.Calc_CPP(1000)

fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
ims = []

#loop over time steps
for t, con in enumerate(C):
    if t % 10 == 0:
        # con is a 3D array at this time step: con[s][y][x]
        im1 = ax1.imshow(con[0], animated=True)  # CO
        im2 = ax2.imshow(con[1], animated=True)  # CO2
        im3 = ax3.imshow(con[2], animated=True)  # O2
        ims.append([im1, im2, im3])
    else: 
        continue
ani = ArtistAnimation(fig, ims, 100)
plt.show()

#/Users/yanni/.pyenv/versions/3.12.2/bin/python -m cProfile -o Animation.prof /Users/yanni/VSC_Projects_Folder/testenv/PPP/ppp/Animation_CPP