# This is an attempt to numerically calculate the magnetic field surounding a
# coil with an oscillating current in it.
# Kipling Crossing

import numpy as np
#import matplotlib.pyplot as plt


# Geometry of coil

pi = np.pi
# The law of Biot-savant - an approximation

permeability = pi*4*10**(-7)
#I = 1       # current in Amps (A)


# Fixed Variables
V = 12      # max voltage in sine wave
f = 10000   # Hz
b = 20      # coil length in mm
radiusT = 50      # Radius in mm
turns = 100 # NUmber of turns in the transmitter coil
n=100       #number of slices

def radius_current():
    # input the radius (a) in mm
    L = ((((radiusT/25.4)**2)*(turns**2))/(9*(radiusT/25.4)+10*(b/25.4)))/1000000   # inductance in H (Henerys)

    XL = 2*np.pi*f*L

    R = (87.85/1000000)*2*np.pi*a*n

    Z = ((R**2)+(XL**2))**0.5
    I = V/Z
    # Returns Current in amps
    return I




def small_flux_dencity(dl_point, dl, focus_point):
    I = radius_current()
    dj = np.subtract(focus_point, dl_point)

    dj_hat = dj/np.linalg.norm(dj)
    dB = np.dot((permeability / (4.0 * pi)) * I*turns/(np.linalg.norm(dj)**2), np.cross(dl, dj_hat))

    return dB


#r = input("Radius: ")
#n = input("Number of slices: ")

#r = float(r)
#n = int(n)



def Flux_dencity(rx,ry,rz,r,turns):
    dl = 2 * r * np.tan(pi / n)
    d = [0, 0, -1]		# controles current direction 
    j = [rx, ry, rz]           # focus point
    
    sum = [0,0,0]
    for i in range(n):
        px = r*np.cos(i*(2*pi/n))
        py = r*np.sin(i*(2*pi/n))
    
        v = [px, py, 0]
        v_hat = v/np.linalg.norm(v)
        dl_hat = np.cross(v_hat, d)
    
    
        dl_v = np.dot(dl,dl_hat)
        #print("Radius vector " + str(v)+" dl: " + str(dl_v))
    
        b = small_flux_dencity(v,dl_v,j,r,turns)
        #print(b)
        sum = np.add(sum,b)

    #print("Sum: " + str(sum))
    return sum

print(Flux_dencity(0,0,0.5,100,1))

f = open('output3.csv', 'w')

for i in range(1,1000):
    radius = i/1000.0       # m
    [q,w,response] = Flux_dencity(0, 0, -0.5, radius,1)
    f.write(str(radius) + "," + str(response)+"\n")
