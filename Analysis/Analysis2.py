import numpy as np 
import matplotlib.pyplot as plt
import os,sys
from Tools import *



	





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
plt.xlim([1.0,3.0])
plt.ylim([-2.5,-1.95])
plt.xlabel("T")
plt.grid(1)
#plt.savefig("F.jpg",format="jpeg")

plt.figure(2)
plt.title("lnZ")
plt.plot(dat[:,2],dat[:,3],"o")


newTs, Es = GetE(dat[:,2],dat[:,3])
plt.figure(3)
plt.title("E")
plt.plot(newTs, Es,"o")
plt.xlim([1.0,3.0])
plt.ylim([-2.05,-0.795])
plt.grid(1)
plt.xlabel("T")
#plt.savefig("E.jpg",format='jpeg')

newTs, Cvs = GetCv(dat[:,2],dat[:,3])
plt.figure(4)
plt.title("Cv")
plt.plot(newTs, Cvs,"o")

newTs, dFs = Derivative(dat[:,2],dat[:,4],2,6)
plt.figure(5)
plt.title("dF")
plt.plot(newTs, dFs,"o")
plt.show()

