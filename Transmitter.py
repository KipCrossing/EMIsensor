# This is an attempt to numerically calculate the magnetic field surounding a
# coil with an oscillating current in it.
# Kipling Crossing


import numpy as np
from vpython import *

# Geometry of coil

pi = np.pi

r = input("Radius: ")
n = input("Number of slices: ")

r = float(r)
n = int(n)

dl = 2*r*np.tan(pi/n)

v = []

for i in range(n):
    px = r*np.cos(i*(2*pi/n))
    py = r*np.sin(i*(2*pi/n))

    v[i] = vector(px,py,0)

    print("px: "+str(px)+" py: "+str(py))
    print(v[i])

print("dl: " + str(dl))
print("actual circumference: " + str(pi*r*2) + "   Approximate circumference: " + str(n*dl))


# The law of Biot-savant - an approximation10

permeability = pi*4*10**(-7)

I = 1       # current in Amps (A)

radius = 0.1        # Radius of Transmitter coil in (m)

B = (permeability/(4.0*pi))*I
print("B = " + str(B))
