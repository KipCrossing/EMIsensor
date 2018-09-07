import numpy as np
from Euler import rotation_matrix
import time
import csv

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


RTx = 0.015  # Tx radius (m)

m = 20
Re = 1.0 / (2 * float(m))
y_offset = 0.0

configurations = [(([0, 0, 1], [1 / m * 2, 1 / m * 2, 1 / m * 2]), ([0, 0, 1], [1 / m * 2, 1 / m * 2, 0.15 + 1 / m * 2]), (
[0, 0, 1], [1 + 1 / m * 2, 1 / m * 2, 1 / m * 2]))]  # A list containing the configurations wanting to be analysed

for config in configurations:
    (Tx1, Tx2, Rx) = config  # Tx is is the transmitter coils, Rx is the receiver coils.
    (Rx_axis, Rx_position) = Rx
    (Axis1, Tx_position1) = Tx1
    Hp = Flux_dencity(Tx1, Tx2, RTx, Rx_position)
    print("tick")



    '''
    f = open('2D-'+str(y_offset)+'_'+str(Tx1)+'_'+str(Tx2)+'_'+str(Rx)+'.csv', 'w')
    f.write("x,y,result\n")


    for k in range(-1 * m, 0):
        for i in range(-1 * m, 2 * m + 1):
            r = [float(i) / m, y_offset, float(k) / m]


            try:
                Hr = Flux_dencity(Tx1, Tx2, RTx, r)


                #Hs = Flux_dencity((Hr, r), None, Re, Rx_position)  # Check for magnitude of Hr in Flux_dencity()
                holder = np.cross(r,Hr)
                dL_hat = holder/np.linalg.norm(holder)
                Current = np.dot(np.linalg.norm(Hr), dL_hat)

                Hs = np.dot(np.linalg.norm(Hr),small_flux(r, Current,Rx_position))
                result = (np.linalg.norm(Hr) * (np.dot(Rx_axis, Hs) / np.dot(Rx_axis, Hp)))   # Change factor

                #print(str(i / float(m)) + "," + str(k / float(m)) + "," + str(result))
                f.write(str(i / float(m)) + "," + str(k / float(m)) + "," + str(result) + "\n")
            except RuntimeWarning:
                print("Error")
                pass
        print(k)
    f.close()

    # Simple 2D


    f = open('Simple-2D'+'.csv', 'w')
    f.write("x,y,result\n")
    for k in range(-1 * m, 0):
        sum = 0
        for i in range(-1 * m, 2 * m + 1):
            r = [float(i) / m, y_offset, float(k) / m]
            result = (1/(np.linalg.norm(r)**2))*(1/(np.linalg.norm(np.subtract([1, 0, 0], r)**2)))



            if float(i) / m < 0 or float(i) / m > 1:
                result = -result

            result = result*np.linalg.norm(np.dot([0,0,1],np.subtract([1, 0, 0], r) / np.linalg.norm(np.subtract([1, 0, 0], r))))


            sum += result
            f.write(str(i / float(m)) + "," + str(k / float(m)) + "," + str(result) + "\n")
        print(sum)
    f.close()

    '''

    print("tock")

'''
    k_file = open('3D_Current_dL' + str(y_offset) + '_' + str(Tx1) + '_' + str(Tx2) + '_' + str(Rx) + '.csv', 'w')

    f = open('1D_TODAY='+str(m)+'_' + str(Tx1) + '_' + str(Tx2) + '_' + str(Rx) + '.csv', 'w')
    f.write("z,result\n")

    for k in range(-2*m, 0):
        layer = 0
        for i in range(-4 * m, 4 * m):
            for j in range(-4 * m, 4 * m):
                r = [float(i) / m, float(j) / m, float(k) / m]
                [x, y, z] = r
                Hr = Flux_dencity(Tx1, Tx2, RTx, r)     # RTx is the Radius of Tx - r is the position of the point in question
                #print(Hr)
                holder = np.cross(r, Hr)
                dL_hat = holder / np.linalg.norm(holder)
                Current = np.dot(np.linalg.norm(Hr),dL_hat)
                #  print(str(r)+" - "+str(np.linalg.norm(Hr)))
                [u, v, w] = Current

                k_file.write("%s ,%s ,%s ,%s ,%s ,%s \n" % (x, y, z, u, v, w))

                Hs = np.dot(np.linalg.norm(Hr),small_flux(r, Current,Rx_position))


                layer += np.linalg.norm(Hr)* (np.dot(Rx_axis, Hs) / np.dot(Rx_axis, Hp))  # Change factor - be sure that we include ohms law in the calculations


        print(str(k / float(m))+ ', ' +str(layer*10**8))
        f.write(str(k / float(m)) + "," + str(layer*10**8) + "\n")
    f.close()
'''

plot_vector = []

with open('3D_Current_dL' + str(y_offset) + '_' + str(Tx1) + '_' + str(Tx2) + '_' + str(Rx) + '.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    for row in readCSV:
        try:
            row_float = [float(i) for i in row]
            plot_vector.append(row_float)
        except:
            pass
    plot_vector.pop(0)


the_z = 0
sum = None
B_sum = [0,0,0]
for [x, y, z, u, v, w] in plot_vector:
    coord = [x, y, z]
    current = [u, v, w] 
    #Rx_position = [1 + 1 / m * 2, 1 / m * 2, 0.5 + 1 / m * 2]
    r = np.subtract(Rx_position, coord)
    unit_r = r / np.linalg.norm(r)
    B = np.cross(current, unit_r) / (np.linalg.norm(r)) ** 2
    B_sum = np.add(B_sum , B)
    #response = np.dot(Rx_axis, B)
    if the_z == z:
        response = np.dot(Rx_axis, B_sum)
        sum = response
    else:
        print("%s, %s" % (the_z, sum))
        response = np.dot(Rx_axis, B_sum)
        sum = 0
        B_sum = [0, 0, 0]
        the_z = z
        sum += response

    #print(response)