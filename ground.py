import numpy as np
from Euler import rotation_matrix
import Transmitter

permeability = np.pi * 4 * 10 ** (-7)

turns = 1

r = 0.001

Emf = np.linalg.norm(Transmitter.Flux_dencity(0.5,0,0.25))*(np.pi*r**2)     # 90 degrees offset!!!!!!!!!!!!!!!!!!!

I = Emf / 0.1  # current in Amps (A)

n = 100

dl = 2 * r * np.tan(np.pi / n)

def small_flux_dencity(dl_point, dl, focus_point):


    dj = np.subtract(focus_point, dl_point)

    dj_hat = dj/np.linalg.norm(dj)
    dB = np.dot((permeability / (4.0 * np.pi)) * I*turns/(np.linalg.norm(dj)**2), np.cross(dl, dj_hat))

    return dB




def start_radius_vector(B,R):
    [a,b,c] = B

    [a, b, c] = [float(a), float(b), float(c)]

    try:
        x = (R ** 2 / (1 + ((a ** 2) / (c ** 2)))) ** 0.5
    except ZeroDivisionError:
        x = 0
    try:
        z = (R ** 2 / (1 + ((b ** 2) / (a ** 2)))) ** 0.5
    except ZeroDivisionError:
        z = 0
    return [x,0,z]




def Flux_dencity(rx, ry, rz):
    d = [0, 0, -1]  # controles current direction
    j = [rx, ry, rz]  # focus point

    sum = [0, 0, 0]
    for i in range(n):
        g = start_radius_vector(d,r)  # Start vector
        Axis = Transmitter.Flux_dencity(0.5,0,0.25)

        v = np.dot(rotation_matrix(Axis, i * (2 * np.pi / n)), g)


        v_hat = v / np.linalg.norm(v)
        dl_hat = np.cross(v_hat, d)

        dl_v = np.dot(dl, dl_hat)


        b = small_flux_dencity(v, dl_v, j)

        sum = np.add(sum, b)
    return sum

print(Flux_dencity(1,0,0))