# This is an attempt to numerically calculate the magnetic field surounding a
# coil with an oscillating current in it.
# Kipling Crossing


import numpy as np


# Geometry of coil

pi = np.pi
# The law of Biot-savant - an approximation

permeability = pi*4*10**(-7)
I = 1       # current in Amps (A)
turns = 1000
def small_flux_dencity(dl_point, dl, focus_point):
    dj = np.subtract(focus_point, dl_point)

    dj_hat = dj/np.linalg.norm(dj)
    dB = np.dot((permeability / (4.0 * pi)) * I*turns/(np.linalg.norm(dj)**2), np.cross(dl, dj_hat))

    return dB


#r = input("Radius: ")
#n = input("Number of slices: ")

#r = float(r)
#n = int(n)

r = 1
n = 100

dl = 2*r*np.tan(pi/n)

print("dl: "+ str(dl))

def Flux_dencity(rx,ry,rz):

    d = [0, 0, -1]		# controles current direction 
    j = [rx, ry, rz]           # focus point
    
    sum = [0,0,0]
    for i in range(n):
        px = r*np.cos(i*(2*pi/n))
        py = r*np.sin(i*(2*pi/n))
    
        v = [px, py, 0]
        v_hat = v/np.linalg.norm(v)
        dl_hat = np.cross(v_hat, d)
    
    
        dl_v = np.dot(dl,dl_hat)
        #print("Radius vector " + str(v)+" dl: " + str(dl_v))
    
        b = small_flux_dencity(v,dl_v,j)
        #print(b)
        sum = np.add(sum,b)

    #print("Sum: " + str(sum))
    return sum

print(Flux_dencity(2,0,1))

#print("dl: " + str(dl))
#print("actual circumference: " + str(pi*r*2) + "   Approximate circumference: " + str(n*dl))


#B = ((permeability*I*turns*(r**2)) /(2.0*((r**2 + 0**2)**0.5)**3))		# where 0 is x 
#print(B)
