import numpy as np
from Euler import rotation_matrix
import time



n = 32            # number of slices

EC = 25
permeability = np.pi * 4 * 10 ** (-7)
frequency = 10000       # Hz
turns = 1

#+print((np.linalg.norm(S)**2)*EC*permeability*frequency*2*np.pi/4)


def small_flux(dl_point, dl, focus_point):

    dj = np.subtract(focus_point, dl_point)

    dj_hat = dj/np.linalg.norm(dj)      # Doesn't work for z = 0!!! So don't consider this case.

    dB = np.dot(1.0/(np.linalg.norm(dj)**2),np.cross(dl, dj_hat))

    return dB



def Flux_dencity(Tx1, Tx2, RTx, r):
    (Axis1, Tx_position1) = Tx1


    dl = 2 * RTx * np.tan(np.pi / n)

    Axis_hat = Axis1 / np.linalg.norm(Axis1)

    c = np.cross(Axis1,[1,0,0])
    g = (c/ np.linalg.norm(c))*RTx
    #print(g)


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


RTx = 0.015  # Tx radius

m = 20
Re = 1.0/(2*float(m))

configuration = [(([0,-1,0],[0,0,RTx]),None,([0,-1,0],[1,0,RTx])),(([0,0,1],[0,0,0]),None,([0,0,1],[1,0,0]))]  #,(([-1,0,1],[0,0,0]),None,([1,0,1],[1,0,0]))]

for tuple in configuration:
    (Tx1,Tx2,Rx) = tuple
    (Rx_axis, Rx_position) = Rx
    Hp = Flux_dencity(Tx1, Tx2, 0.05, Rx_position)
    '''
    f = open('2D-0.2_'+str(Tx1)+'_'+str(Tx2)+'_'+str(Rx)+'.csv', 'w')
    f.write("x,y,result\n")

    for k in range(-2 * m, 0):
        for i in range(-1 * m, 2 * m + 1):
            r = [float(i) / m, -0.2, float(k) / m]
            try:
                Hr = Flux_dencity(Tx1, Tx2, 0.05, r)
                Hs = Flux_dencity((Hr, r), None, Re, Rx_position)  # Check for magnitude of Hr in Flux_dencity()

                result = (np.linalg.norm(Hr) * (np.dot(Rx_axis, Hs)*Re**2 / np.dot(Rx_axis, Hp)))   # Change factor
                # print(str(i / float(m)) + "," + str(k / float(m)) + "," + str(result))
                f.write(str(i / float(m)) + "," + str(k / float(m)) + "," + str(result) + "\n")
            except RuntimeWarning:
                print("Error")
                pass
        print(k)
    f.close()

    '''


    f = open('1D_0.1_' + str(Tx1) + '_' + str(Tx2) + '_' + str(Rx) + '.csv', 'w')
    f.write("x,y,result\n")

    for k in range(-1 * m, 0):
        layer = 0
        for i in range(-1 * m, 2 * m + 1):
            for j in range(-1 * m, 1 * m+1):
                r = [float(i) / m, float(j) / m, float(k) / m]
                try:
                    Hr = Flux_dencity(Tx1, Tx2, 0.05, r)
                    Hs = Flux_dencity((Hr, r), None, Re, Rx_position)  # Check for magnitude of Hr in Flux_dencity()

                    layer += np.linalg.norm(Hr)* (np.dot(Rx_axis, Hs)*Re**2 / np.dot(Rx_axis, Hp))  # Change factor
                    # print(str(i / float(m)) + "," + str(k / float(m)) + "," + str(result))

                except RuntimeWarning:
                    print("Error")
                    pass
        print(str(k / float(m))+ ', ' +str(layer))
        f.write(str(k / float(m)) + "," + str(layer) + "\n")
    f.close()
