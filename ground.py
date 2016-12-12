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



def Flux_dencity_ground(rx, ry, rz):
    # Inputs the position vector at where the reading is being focused
    receiver_position = [1, 0, 0]  # focus point

    Flux = [0, 0, 0]

    Axis = Transmitter.Flux_dencity(rx, ry, rz)

    Axis_hat = Axis/np.linalg.norm(Axis)

    ground_position = [rx, ry, rz]

    g = start_radius_vector(Axis, r)  # Start vector

    dl = 2 * r * np.tan(np.pi / n)

    Emf = np.linalg.norm(Axis) * (np.pi * r ** 2)  # 90 degrees offset!!!!! (in time)

    I = Emf / 0.1  # current in Amps (A)

    for i in range(n):

        v = np.dot(rotation_matrix(Axis, i * (2 * np.pi / n)), g)

        dl_point = np.add(ground_position, v)

        v_hat = v / np.linalg.norm(v)
        dl_hat = np.cross(v_hat, Axis_hat)

        dl_v = np.dot(dl, dl_hat)

        b = small_flux_dencity(dl_point, dl_v, receiver_position, I)

        Flux = np.add(Flux, b)
    return Flux

'''
f = open('output_y0.0.csv', 'w')
f.write("x,y,vertical,horizontal\n")

for i in range(-50, 101):
    for k in range(-100, -10):
        result = Flux_dencity_ground(i/50.0,0.0,k/50.0)
        if result[2]*(10**20) > -10.0 and result[2]*(10**20) < 10.0:
            f.write(str(i / 50) + "," + str(k / 50) + "," + str(result[0]*(10**20)) + "," + str(result[2]*(10**20))+ "\n")

    print(i/50.0)

f.close()
'''



fi = open('output_later.csv', 'w')

fi.write("z,layer\n")

for k in range(-30, 1):
    layer = 0
    for i in range(-10, 21):
        for j in range(-10,10):
            result = Flux_dencity_ground(i/10.0,j/10.0,k/10.0)
            layer += result[2]
    fi.write(str(k / 10) + "," + str(layer) + "\n")
    print(k)
fi.close()
