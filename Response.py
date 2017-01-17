import numpy as np
from Euler import rotation_matrix

Tx = [0,0,1]
Rx = [1,0,0]


EC = 25
permeability = np.pi * 4 * 10 ** (-7)
frequency = 10000       # Hz
S = [1, 0, 0]
turns = 1


RTx = 0.05
n = 1000         # number of slices

r = [0.5, 0, 0.5]

#print((np.linalg.norm(S)**2)*EC*permeability*frequency*2*np.pi/4)


def small_flux(dl_point, dl, focus_point):

    dj = np.subtract(focus_point, dl_point)

    dj_hat = dj/np.linalg.norm(dj)      # Doesn't work for z = 0!!! So don't consider this case.

    dB = np.dot(1.0/(np.linalg.norm(dj)**2),np.cross(dl, dj_hat))

    return dB



def Flux_dencity(Axis, Tx_position, RTx, r):
    (Axis1, Axis2) = Axis
    (Tx_position1, Tx_position2) = Tx_position

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

    if Axis2 != None:

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




f = open('OUTPUT.csv', 'w')
f.write("x,y,result\n")

m = 10
Rx_position = [1,0,0]
Re = 1.0/(2*float(m))
Hp = Flux_dencity((Tx, [0,0,-1]),([0,0,0],[2,0,0]),0.05,Rx_position)
bottom = np.dot(Rx,Hp)
for k in range(-2 * m, 0):
    for i in range(-1 * m, 3 * m + 1):
        r = [i,0.0,k]
        Hr = Flux_dencity((Tx, [0,0,-1]),([0,0,0],[2,0,0]),0.05,r)
        Hs = Flux_dencity((Hr, None),(r,None),Re,Rx_position)*10**8     # Check for magnitude of Hr in Flux_dencity()

        result = np.dot(Rx,Hs) #/ np.dot(Rx,Hp)          # np.linalg.norm(Hr)*
        #print(str(i / float(m)) + "," + str(k / float(m)) + "," + str(Hr))
        f.write(str(i / float(m)) + "," + str(k / float(m)) + "," + str(result) + "\n")

f.close()

