import numpy as np
from Euler import rotation_matrix
import time
import csv
import pandas as pd

n = 32  # number of slices

EC = 25
permeability = np.pi * 4 * 10 ** (-7)
frequency = 10000  # Hz
turns = 1


# +print((np.linalg.norm(S)**2)*EC*permeability*frequency*2*np.pi/4)


def small_flux(dl_point, dl, focus_point):
    dj = np.subtract(focus_point, dl_point)

    dj_hat = dj / np.linalg.norm(dj)  # Doesn't work for z = 0!!! So don't consider this case.

    dB = np.dot(1.0 / (np.linalg.norm(dj) ** 2), np.cross(dl, dj_hat))

    return dB


def Flux_dencity(Tx1, Tx2, RTx, r):
    (Axis1, Tx_position1) = Tx1

    dl = 2 * RTx * np.tan(np.pi / n)

    Axis_hat = Axis1 / np.linalg.norm(Axis1)

    c = np.cross(Axis1, [1, 0, 0])

    g = (c / np.linalg.norm(c)) * RTx

    sum = [0, 0, 0]
    for i in range(n):
        v = np.dot(rotation_matrix(Axis1, i * (2 * np.pi / n)), g)

        dl_point = np.add(Tx_position1, v)

        v_hat = v / np.linalg.norm(v)
        dl_hat = np.cross(v_hat, Axis_hat)

        dl_v = np.dot(dl, dl_hat)

        b = small_flux(dl_point, dl_v, r)

        sum = np.add(sum, b)

    if Tx2 != None:
        (Axis2, Tx_position2) = Tx2
        Axis_hat2 = Axis2 / np.linalg.norm(Axis2)

        for i in range(n):
            v = np.dot(rotation_matrix(Axis2, i * (2 * np.pi / n)), g)

            dl_point = np.add(Tx_position2, v)

            v_hat = v / np.linalg.norm(v)
            dl_hat = np.cross(v_hat, Axis_hat2)

            dl_v = np.dot(dl, dl_hat)

            b = small_flux(dl_point, dl_v, r)

            sum = np.add(sum, b)

    return sum


RTx = 0.005  # Tx radius (m)
m = 10000
Re = 1.0 / (2 * float(m))
y_offset = 0.0

configurations = [(([0, 0, 1], [0, 0, 0]), None ,([0, 0, 1], [1,1, 1]))]  # A list containing the configurations wanting to be analysed

for config in configurations:
    (Tx1, Tx2, Rx) = config  # Tx is is the transmitter coils, Rx is the receiver coils.
    (Rx_axis, Rx_position) = Rx
    (Axis1, Tx_position1) = Tx1
    Hp = Flux_dencity(Tx1, Tx2, RTx, Rx_position)
    print(Hp)

x_n = -3
x_p = 0
y_n = -1
y_p = 1
z_n = -2
z_p = 0

Ey = 0
for i in range(x_n*m,x_p*m+1):
    mag = Flux_dencity(Tx1, Tx2, RTx, [float(i)/m-0.00005,0,0])
    Ey += mag[2]/m
    print('X:',float(i)/m,'  -  ',Ey,'  -  ',mag[2]*1/m)

for j in range(y_n*m,y_p*m+1):
    print('Y:',float(j)/m)

for k in range(z_n*m,z_p*m+1):
    print('Z:',float(k)/m)
