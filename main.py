import numpy as np
from Euler import rotation_matrix
import time
#import csv
#import pandas as pd
#%matplotlib inline
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

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


# initial perameters
m = 1000
m_step = 1.0/(m*2)

RTx = 0.005  # Tx radius (m)

Re = 1.0 / (2 * float(m))
y_offset = 0.0

configurations = [(([0, 0, 1], [0, 0, 0]), None ,([0, 0, 1], [1,1,1]))]  # A list containing the configurations wanting to be analysed

for config in configurations:
    (Tx1, Tx2, Rx) = config  # Tx is is the transmitter coils, Rx is the receiver coils.
    (Rx_axis, Rx_position) = Rx
    (Axis1, Tx_position1) = Tx1
    Hp = Flux_dencity(Tx1, Tx2, RTx, Rx_position)
    print(Hp)

'''

x_n = 0
x_p = 1
y_n = 0
y_p = 0
z_n = -1
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

    print('Z:',float(k)/m, 'B', sum)
    z.append([float(k)/m,sum])

z.reverse()
z

#x = []

x_temp = []
y_temp = []
z_temp = []

number = 0

for k in z:
    x_temp = [0]
    z_temp.append([k[1]])
    print(k[0],k[1])
    y_temp.append(k[0])
    #x.append([0.0,k[0],k[1]])
    Electric_vector = k[1]
    for i in range(int(x_n*m),int(x_p*m)):
        r = [float(i)/m+m_step,0,k[0]]
        [float(i)/m+m_step,0,number*0.0001]
        mag = Flux_dencity(Tx1, Tx2, RTx, r)/m
        #print(mag)
        Electric_vector -= mag[2]
        #x.append([float(i)/m+(m_step*2),k[0],Electric_vector])
        x_temp.append(float(i)/m+(m_step*2))
        z_temp[len(z_temp)-1].append(Electric_vector)


len(x_temp)
len(y_temp)
len(z_temp)
z_temp
x



# x_data = []
# y_data = []
# z_data = []
#
#
# for lst in x:
#
#     x_data.append(lst[0])
#     y_data.append(lst[1])
#     z_data.append(lst[2])


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
X, Y = np.meshgrid(x_temp, y_temp)
Z = np.array(z_temp)


ax.plot_surface(X, Y, Z, alpha=0.5, rstride=1, cstride=1)

plt.show()





ax = plt.axes(projection='3d')
ax.scatter(x_data, z_data, z_data)  #,cmap='viridis', edgecolor='none');


# import scipy.interpolate as interp
#
# plotx,ploty, = np.meshgrid(np.linspace(np.min(X),np.max(X),10),\
#                            np.linspace(np.min(Y),np.max(Y),10))
# plotz = interp.griddata((X,Y),Z,(plotx,ploty),method='linear')
#
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(plotx,ploty,plotz,cstride=1,rstride=1,cmap='viridis')  # or 'hot'




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

#z.reverse()

Bz
len(Bz)
len(x)
len(y)

z.reverse()
y.pop(1)
y.pop(500)
z
y.append(0)
y.insert(0,z[0])
len(z)

Bz_check = []
for j in range(0,len(y)-1):
    print(j)
    print(y[j+1] - y[j])
    Bz_check.append(-y[j+1] + y[j])
#y.pop(0)
len(Bz_check)
y

Bz_check.append(0)

for a in range(len(Bz_check)):
    print(round(Bz_check[a],8) == round(Bz[a],8))

#Adding axis labels
plt.xlabel('x axis')
plt.ylabel('Electric field')
#adding title
plt.title('Ey along the X-asis')


plt.plot(x,y,'-',c='red',label = "Ey on x")
plt.plot(x,z,'-',c='purple',label = "Exy on z")
plt.plot(x,Bz,'-',c='blue',label = "dEy/dx on x")
#plt.plot(x,Bz_check,'-',c='blue',label = "dEy/dx on x check")
plt.plot(0.005,0,'o',label = "coil into the page")
plt.plot(-0.005,0,'o',label = "coil out of the page")
plt.legend(loc ='upper right')
plt.show()

# for j in range(y_n*m,y_p*m+1):
#     print('Y:',float(j)/m)


sum = 0
for k in range(int(z_n*m),int(z_p*m+1)):
    sum += Flux_dencity(Tx1, Tx2, RTx, [0,0,float(k)/m-m_step])[2]/(m)
    print('Z:',float(k)/m, 'B', sum)

dx = 0.0001
dz = 0.0001

y[1]
z[0]

#Flux_dencity(Tx1, Tx2, RTx, [m_step,0,-m_step])/(m)
Flux_dencity(Tx1, Tx2, RTx, [0.0001,0,-m_step])/(m)
Flux_dencity(Tx1, Tx2, RTx, [0,0,0.0001])/(m)
print("start")

number = 0
E1 = z[number]
for i in range(int(x_n*m),int(0.1*m)):
    mag = Flux_dencity(Tx1, Tx2, RTx, [float(i)/m+m_step,0,number*0.0001])/(m)     #0.0001
    E1 -= mag[2]
    print('X:',float(i+1)/m,'  -  ',E1,'  -  ',mag[2]*1/m)

sum = 0
for k in range(int(z_n*m),int(z_p*m+1)):
    #sum +=Flux_dencity(Tx1, Tx2, RTx, [0,0,float(k)/m-m_step])[2]/(m)
    sum += Flux_dencity(Tx1, Tx2, RTx, [0.0001,0,float(k)/m-m_step])[2]/(m)
    print("Z:",k,float(k)/m,"Ey:",sum)


#z[1]-Flux_dencity(Tx1, Tx2, RTx, [m_step,0,0.0001])[2]/(m)
'''
