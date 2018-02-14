import numpy as np 
import matplotlib.pyplot as plt
import os,sys

def Div(x,y):
    newx = (x[1:] + x[0:-1])/2

    dydx = ( y[1:] - y[0:-1] ) / (x[1:] - x[0:-1]) 
    
    return dydx, newx


def GetE(lnZs,Ts):
    dlnZdT, newT = Div(Ts,lnZs) 
    return newT**2 * dlnZdT , newT


def GetCv(lnZs,Ts):
    
    dlnZdT  , newT = Div(Ts,lnZs)
    dlnZdT2 , newT = Div(newT, dlnZdT)
 
    return newT**2 * dlnZdT2 , newT





f = open(sys.argv[1],"r")
lines = f.readlines()
f.close()
dat = []
for line in lines:
    tmp = np.array(line.strip().split(" "),dtype=np.float)
    dat.append(tmp)
dat = np.array(dat)


plt.figure(1)
plt.title("F")
plt.plot(dat[:,2],dat[:,4],"o")

plt.figure(2)
plt.title("lnZ")
plt.plot(dat[:,2],dat[:,3],"o")


Es , newTs = GetE(dat[:,3],dat[:,2])
plt.figure(3)
plt.title("E")
plt.plot(newTs, Es,"o")

Cvs , newTs = GetCv(dat[:,3],dat[:,2])
plt.figure(4)
plt.title("Cv")
plt.plot(newTs, Cvs,"o")
plt.show()
