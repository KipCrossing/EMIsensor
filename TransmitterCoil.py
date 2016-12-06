import numpy as np
import matplotlib.pyplot as plt


# Fixed Variables
V = 12      # max voltage in sine wave
f = 0   # Hz
b = 20      # coil length in mm


def radius_current(a,n):
    # input the radius (a) in mm
    L = ((((a/25.4)**2)*(n**2))/(9*(a/25.4)+10*(b/25.4)))/1000000   # inductance in H (Henerys)

    XL = 2*np.pi*f*L

    R = (87.85/1000000)*2*np.pi*a*n

    Z = ((R**2)+(XL**2))**0.5
    I = V/Z
    # Returns Current in amps
    return I



