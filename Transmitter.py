# This is an attempt to numerically calculate the magnetic field surounding a
# coil with an oscillating current in it.
# Kipling Crossing

import numpy as np
#import matplotlib.pyplot as plt
from Euler import rotation_matrix


# Geometry of coil

pi = np.pi
# The law of Biot-savant - an approximation

permeability = pi*4*10**(-7)
#I = 1       # current in Amps (A)


# Fixed Variables
V = 12      # max voltage in sine wave
f = 10000   # Hz
radiusT =0.45      # Radius in m
LengthT = radiusT   # coil length in mm (L >= 0.8r)
turns = 500 # NUmber of turns in the transmitter coil
n = 1000       #number of slices

def radius_current():
    # input the radius (a) in mm
    L = ((((radiusT/0.0254)**2)*(turns**2))/(9*(radiusT/0.0254)+10*(LengthT/0.0254)))/1000000   # inductance in H (Henerys)

    XL = 2*np.pi*f*L

    R = (87.85/1000)*2*np.pi*radiusT*n

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




def Flux_dencity(rx,ry,rz, tuple):

    (Axis, Axis2) = tuple

    dl = 2 * radiusT * np.tan(pi / n)

    Axis_hat = Axis/np.linalg.norm(Axis)
    Transmitter_position = [0,0,0]

    g = [0,radiusT,0]

    j = [rx, ry, rz]           # focus point
    
    sum = [0,0,0]
    for i in range(n):
        v = np.dot(rotation_matrix(Axis, i * (2 * np.pi / n)), g)

        dl_point = np.add(Transmitter_position, v)

        v_hat = v / np.linalg.norm(v)
        dl_hat = np.cross(v_hat, Axis_hat)

        dl_v = np.dot(dl, dl_hat)

        b = small_flux_dencity(dl_point, dl_v, j)

        sum = np.add(sum, b)




    if Axis2 != None:

        Axis_hat2 = Axis2 / np.linalg.norm(Axis2)
        Transmitter_position2 = [2, 0, 0]

        for i in range(n):
            v = np.dot(rotation_matrix(Axis2, i * (2 * np.pi / n)), g)

            dl_point = np.add(Transmitter_position2, v)

            v_hat = v / np.linalg.norm(v)
            dl_hat = np.cross(v_hat, Axis_hat2)

            dl_v = np.dot(dl, dl_hat)

            b = small_flux_dencity(dl_point, dl_v, j)

            sum = np.add(sum, b)

    return sum





'''
file = open('contureplot.csv', 'w')
file.write("x,y,response\n")
for i in range(-10,11):
    for k in range(1,11):
        result = Flux_dencity(i/10,0,k/10.0)

        print(result)
        #print(str(i/10)+","+str(k/10)+","+str(result[0]) + ","+str(result[2]))
        #file.write(str(i/10)+","+str(k/10)+","+str(np.linalg.norm(result)) + "\n")

file.close()
'''

'''
f = open('output3.csv', 'w')

for i in range(1,1000):
      # m
    [q,w,response] = Flux_dencity(0, 0, -0.5, radius,1)
    f.write(str(radius) + "," + str(response)+"\n")
'''
