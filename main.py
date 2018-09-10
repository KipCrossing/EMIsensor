import numpy as np
from Euler import rotation_matrix
import time
import csv
import pandas as pd
import matplotlib.pyplot as plt

n = 128  # number of slices

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



m = 1000
m_step = 1.0/(m*2)
m_step
RTx = 0.005+m_step  # Tx radius (m)
RTx
Re = 1.0 / (2 * float(m))
y_offset = 0.0

configurations = [(([0, 0, 1], [0, 0, 0]), None ,([0, 0, 1], [1,1, 1]))]  # A list containing the configurations wanting to be analysed

for config in configurations:
    (Tx1, Tx2, Rx) = config  # Tx is is the transmitter coils, Rx is the receiver coils.
    (Rx_axis, Rx_position) = Rx
    (Axis1, Tx_position1) = Tx1
    Hp = Flux_dencity(Tx1, Tx2, RTx, Rx_position)
    print(Hp)

x_n = -10
x_p = 0.1
y_n = -1
y_p = 1
z_n = -0.01
z_p = 0

dE = [0,0,0]
E = 0
x=[]
y=[]
for i in range(int(x_n*m)+1,int(x_p*m)+1):
    r = [float(i)/m,0,0]
    mag = Flux_dencity(Tx1, Tx2, RTx, r)
    r_hat =r / np.linalg.norm(r)


    dE = np.cross(r_hat,mag)

    E+=dE[1]/m
    print('X:',float(i)/m,'  -  ',E,'  -  ',mag[2]*1/m)
    if -0.05*m<i<0.05*m:
        x.append(float(i)/m)
        y.append(E)


print(Flux_dencity(Tx1, Tx2, RTx, [0,0,0])[2]/(m))

(-6.14586631254/-1.25688945198)*np.pi

-6.14586631254*np.pi
(-1.25688945198*m/np.pi*2)

#Adding axis labels
plt.xlabel('x')
plt.ylabel('r_hat cross Ey')
#adding title
plt.title('Ey along the X-asis')

plt.plot(x,y,'-',c='red')
plt.plot(0.005,0,'o')
plt.plot(-0.005,0,'o')
plt.show()
# for j in range(y_n*m,y_p*m+1):
#     print('Y:',float(j)/m)
#
# for k in range(int(z_n*m),z_p*m+1):
#     print('Z:',float(k)/m, 'B', Flux_dencity(Tx1, Tx2, RTx, [0,0,float(k)/m])/(m))
