import numpy as np 
import matplotlib.pyplot as plt
import os,sys
from Tools import *


## Read File:
Mdir = sys.argv[1]
if not os.path.exists(Mdir):
    print("%s not exists."%(Mdir))
    exit(1)


## list files:
sname = [x for x in os.listdir(Mdir) if not os.path.isdir(os.path.join(Mdir,x)) and not 'result' in x]



rawR = []
for fname in sname:
    print(fname)
    #if not '256' in fname:
    #    continue
    f = open(os.path.join(Mdir,fname),"r")
    lines = f.readlines()
    f.close()
    dat = []
    for line in lines:
        tmp = np.array(line.strip().split(" "),dtype=np.float)
        dat.append(tmp)
    dat = np.array(dat)
    #rawR.append(dat)

    newTs,Es = GetE_simple(dat[:,2],dat[:,3]) 
    f = open(os.path.join(Mdir,fname+'.resultE'),'w')
    for a in range(len(newTs)):
        f.write("%11.11lf %11.11lf %11.11lf %11.11lf\n"%(dat[a,0],dat[a,1],newTs[a],Es[a]))
    f.close()
    


    ## Save:
    newTs,Cvs = GetCv_simple(dat[:,2],dat[:,3])
    f = open(os.path.join(Mdir,fname+'.resultCv'),'w')
    for a in range(len(newTs)):
        f.write("%11.11lf %11.11lf %11.11lf %11.11lf\n"%(dat[a,0],dat[a,1],newTs[a],Cvs[a]))
    f.close()

#rawR = np.array(rawR)
 

"""
plt.figure(1)
plt.title("F")
plt.plot(dat[:,2],dat[:,4],"o")
plt.xlabel("T")
plt.ylabel("F")
plt.grid()
plt.savefig("F.jpg",format='jpeg',dpi=200)

plt.figure(2)
plt.title("lnZ")
plt.plot(dat[:,2],dat[:,3],"o")

Es , newTs = GetE(dat[:,3],dat[:,2])
plt.figure(3)
plt.title("E")
plt.plot(newTs, Es,"o")
plt.xlabel("T")
plt.ylabel("E")
plt.grid()
plt.savefig("E.jpg",format='jpeg',dpi=200)

Cvs , newTs = GetCv(dat[:,3],dat[:,2])
plt.figure(4)
plt.title("Cv")
plt.plot(newTs, Cvs,"o")
plt.xlabel("T")
plt.ylabel("Cv")
plt.savefig("Cv.jpg",format='jpeg',dpi=200)
plt.grid()

plt.show()
"""
