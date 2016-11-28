import math
import numpy as np
import matplotlib.pyplot as plt


def Emf(n, area, dH, dt):
    mu0 = 4*math.pi*10**(-7)

    v = -mu0*n*area*(dH/dt)
    return v

D = 0.1     # m
Di = 0.09   # m

A = (math.pi/8)*(D+Di)**2
print(A)
print(Emf(10000, A, 20, 0.5))




X = np.linspace(0,1,360,endpoint=True)

Bm = 100
frequency = 10
om = 2*np.pi*frequency
Y = Bm*np.sin(om*X)
plt.plot(X, Y)
plt.xlabel('time (s)')
plt.ylabel('Emf (volts)')
plt.show()
