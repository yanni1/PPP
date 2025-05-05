import Conc_Calc_CPP
import random
import matplotlib.pyplot as plt
import params_module


vp = Conc_Calc_CPP.CC_CPP(100)
#print(vp)


p = params_module.params(100)

#print(p.v)


plt.figure(figsize=(8, 5))
plt.plot(p.x, p.v, label="Velocity Profile", color="blue", linewidth=2)
plt.xlabel("x Position")
plt.ylabel("Velocity")
plt.title("Velocity Profile")
plt.grid(True)

# Show the plot
plt.show()