import numpy as np
from Euler import rotation_matrix
import Transmitter

permeability = np.pi * 4 * 10 ** (-7)

turns = 1

r = 0.01       # radius of eddy currents in m

n = 100         # number of slices


def small_flux_dencity(dl_point, dl, focus_point, I):

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
    # Inputs the position vector at where the reading is being focused
    receiver_position = [1, 0, 0]  # focus point

    Flux = [0, 0, 0]

    Axis = Transmitter.Flux_dencity(rx, ry, rz)

    g = start_radius_vector(Axis, r)  # Start vector

    dl = 2 * r * np.tan(np.pi / n)

    Emf = np.linalg.norm(Axis) * (np.pi * r ** 2)  # 90 degrees offset!!!!! (in time)

    I = Emf / 0.1  # current in Amps (A)

    for i in range(n):

        v = np.dot(rotation_matrix(Axis, i * (2 * np.pi / n)), g)

        dl_point = np.add(Axis, v)

        v_hat = v / np.linalg.norm(v)
        dl_hat = np.cross(v_hat, Axis)

        dl_v = np.dot(dl, dl_hat)

        b = small_flux_dencity(dl_point, dl_v, receiver_position, I)

        Flux = np.add(Flux, b)
    return Flux
for i in range(2,10):
    for k in range(10):
        result = Flux_dencity(i/10,0,k/10.0)
        print(result[2])




'''
f = open('output.csv', 'w')
j=0
for i in range(0,10):
    print(i)
    for k in range(0, 10):
        result = Flux_dencity(i/10.0,j/10.0,k/10.0)
        #print(str(i) + "," +str(j) + "," +str(k))
        f.write(str(i/10.0) +  "," +str(k/10.0) + ","  + str(result[2]) + "\n")
        #print(str(i/10.0) + "," +str(j/10.0) + "," +str(k/10.0) + "," +str(result[0]) + "," +str(result[1]) + "," +str(result[2]) + "\n")
f.close()
'''