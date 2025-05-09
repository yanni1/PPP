# Animation update function
from calc import calc
import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation
#import Conc_Calc_CPP

C = calc(25)
fig, (ax1, ax2,ax3) = plt.subplots(1,3)
ims=[]
frame = 0
for con in C:
    frame=frame+1
    if frame % 10 == 0:
        im1 = ax1.imshow(con[0], animated=True)
        im2 = ax2.imshow(con[1], animated=True)
        im3 = ax3.imshow(con[2], animated=True)
        ims.append([im1,im2,im3])
    else:
        continue
ani = ArtistAnimation(fig, ims, 100)
plt.show()