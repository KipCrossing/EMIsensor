import numpy as np
from Euler import rotation_matrix
import Transmitter
import _thread
import time

permeability = np.pi * 4 * 10 ** (-7)

turns = 1

r = 0.01       # radius of eddy currents in m
n = 1000         # number of slices


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



def Flux_dencity_ground(rx, ry, rz, tuple):
    # Inputs the position vector at where the reading is being focused
    receiver_position = [1, 0, 0]  # focus point

    Flux = [0, 0, 0]

    Axis = Transmitter.Flux_dencity(rx, ry, rz, tuple)

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





configuration  = [([0,0,1],[0,0,1]),([0,0,1],[0,0,-1]),([0,0,1],None),([0,1,0],[0,1,0]),([0,1,0],[0,-1,0]),([0,-1,0],None),([1,0,0],[1,0,0]),([1,0,0],[-1,0,0]),([1,0,0],None),([1, 0, 1], [-1, 0, 1])]


for tuple in configuration:

    start = time.time()
    print(str(tuple))


    f = open('1Ddepth/'+str(tuple)+'.csv', 'w')
    for k in range(-10, 1):
        layerx = 0
        layery = 0
        layerz = 0
        for i in range(-10, 31):
            for j in range(-2, 3):
                result = Flux_dencity_ground(i / 10.0, j / 10.0, k / 20.0,tuple)
                layerx += result[0]
                layery += result[1]
                layerz += result[2]
        f.write(str(k / 20.0) + "," + str(layerx)+ "," + str(layery)+ "," + str(layerz) + "\n")
    f.close()


    f = open('2Dfield/0m/' + str(tuple) + '.csv', 'w')
    f.write("x,y,vertical,horizontal\n")
    m = 10
    for i in range(-1 * m, 3 * m + 1):
        for k in range(-2 * m, 0):
            result = Transmitter.Flux_dencity(i / float(m), 0.0, k / float(m), tuple)
            f.write(str(i / float(m)) + "," + str(k / float(m)) + "," + str(result[0]) + "," + str(
                result[2] ) + "\n")
    f.close()



    f = open('2Dfield/0.2m/' + str(tuple) + '.csv', 'w')
    f.write("x,y,vertical,horizontal\n")
    m = 10
    for i in range(-1 * m, 3 * m + 1):
        for k in range(-2 * m, 0):
            result = Transmitter.Flux_dencity(i / float(m), 0.2, k / float(m),tuple)
            f.write(str(i / float(m)) + "," + str(k / float(m)) + "," + str(result[0] ) + "," + str(
                result[2] ) + "\n")
    f.close()



    f = open('2Dresponse/0m/' + str(tuple) + '.csv', 'w')

    f.write("x,y,Rx,Ry,Rz\n")
    m = 20
    for i in range(-1 * m, 3 * m + 1):
        for k in range(-2 * m, 1):
            result = Flux_dencity_ground(i / float(m), 0.0, k / float(m),tuple)
            f.write(str(i / float(m)) + "," + str(k / float(m)) + "," + str(result[0] )+ "," + str(result[1] ) + "," + str(
                result[2] ) + "\n")
    f.close()

    f = open('2Dresponse/0.2m/' + str(tuple) + '.csv', 'w')
    f.write("x,y,Rx,Ry,Rz\n")
    m = 20
    for i in range(-1 * m, 3 * m + 1):
        for k in range(-2 * m, 1):
            result = Flux_dencity_ground(i / float(m), 0.2, k / float(m),tuple)
            f.write(str(i / float(m)) + "," + str(k / float(m)) + "," + str(result[0])+ "," + str(result[1] ) + "," + str(
                result[2]) + "\n")
    f.close()

    end = time.time()

    print("Time: "+str((end - start)/60.0)+" minutes")


