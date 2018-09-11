import numpy as np
from Euler import rotation_matrix
import time
import csv
import pandas as pd
import matplotlib.pyplot as plt

n = 1024  # number of slices

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


# initial perameters
m = 10000
m_step = 1.0/(m*2)

RTx = 0.005  # Tx radius (m)

Re = 1.0 / (2 * float(m))
y_offset = 0.0

# configurations = [(([0, 0, 1], [0, 0, 0]), None ,([0, 0, 1], [1,1, 1]))]  # A list containing the configurations wanting to be analysed
#
# for config in configurations:
#     (Tx1, Tx2, Rx) = config  # Tx is is the transmitter coils, Rx is the receiver coils.
#     (Rx_axis, Rx_position) = Rx
#     (Axis1, Tx_position1) = Tx1
#     Hp = Flux_dencity(Tx1, Tx2, RTx, Rx_position)
#     print(Hp)

x_n = 0
x_p = 2
y_n = -1
y_p = 1
z_n = -4
z_p = 0

dE = [0,0,0]
E = 0
x=[]
y=[0]
z=[]
Bz=[]
sum = 0

for k in range(int(z_n*m),int(z_p*m+1)):
    sum +=Flux_dencity(Tx1, Tx2, RTx, [0,0,float(k)/m-m_step])[2]/(m)
    if k%0.1*m == 0:
        print('Z:',float(k)/m, 'B', sum)
    if k > -0.05*m :
        z.append(sum)

start = z[len(z)-1]
E=start
for i in range(int(x_n*m),int(x_p*m)):
    r = [float(i)/m+m_step,0,0]
    mag = Flux_dencity(Tx1, Tx2, RTx, r)
    r_hat =r / np.linalg.norm(r)

    dE = np.cross(r_hat,mag)

    E+=dE[1]/m
    k=i

    #print('Z:',float(k)/m, 'B', sum)


    print('X:',float(i+1)/m,'  -  ',E,'  -  ',mag[2]*1/m)
    if -0.05*m<i<0.05*m:
        x.append(float(i)/m)
        y.append(E)
        Bz.append(mag[2]*1/m)


print('Finished.')
#print(Flux_dencity(Tx1, Tx2, RTx, [0,0,0])[2]/(m))


# Bz
# len(Bz)
#
# len(x)
# len(y)
#
# Bz_check = []
# for j in range(0,len(y)-1):
#     print(j)
#     print(y[j+1] - y[j])
#     Bz_check.append(y[j+1] - y[j])
# y.pop(0)
# len(Bz_check)
#
#
# for a in range(len(Bz_check)):
#     print(round(Bz_check[a],8) == round(Bz[a],8))
#
# #Adding axis labels
# plt.xlabel('x or z axis')
# plt.ylabel('Electric field')
# #adding title
# plt.title('Ey along the X-asis')
#
#
# plt.plot(x,y,'-',c='red',label = "Ey on x")
# plt.plot(x,z,'-',c='purple',label = "Exy on z")
# plt.plot(x,Bz,'-',c='yellow',label = "dEy/dx on x")
# plt.plot(x,Bz_check,'-',c='blue',label = "dEy/dx on x check")
# plt.plot(0.005,0,'o')
# plt.plot(-0.005,0,'o')
# plt.legend(loc ='lower left')
# plt.show()
# # for j in range(y_n*m,y_p*m+1):
# #     print('Y:',float(j)/m)
# #
# sum = 0
# for k in range(int(z_n*m),int(z_p*m+1)):
#     sum +=Flux_dencity(Tx1, Tx2, RTx, [0,0,float(k)/m-m_step])[2]/(m)
#     print('Z:',float(k)/m, 'B', sum)
