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
sname = [x for x in os.listdir(Mdir) if not os.path.isdir(os.path.join(Mdir,x))]



rawR = []
for fname in sname:
    print(fname)
    if not '256' in fname:
        continue
    f = open(os.path.join(Mdir,fname),"r")
    lines = f.readlines()
    f.close()
    dat = []
    for line in lines:
        tmp = np.array(line.strip().split(" "),dtype=np.float)
        dat.append(tmp)
    dat = np.array(dat)
    #rawR.append(dat)
    
    J1set = np.sort(np.unique(dat[:,0]))
    J2set = np.sort(np.unique(dat[:,1]))
    Tset  = np.sort(np.unique(dat[:,2]))
    hset  = np.sort(np.unique(dat[:,3]))
    dat.view('f8,f8,f8,f8,f8,f8').sort(order=['f0','f1','f3','f2'],axis=0)

    # Dim = [NJ1,NJ2,Nh,NT,6]
    dat = dat.reshape((len(J1set)*len(J2set),len(hset),len(Tset),6))


    # get h=0 sector for E:
    dath0 = dat[:,np.argwhere(hset==0)[0,0],:,:]
        
    datM0 = []
    for j1j2 in range(len(J1set)*len(J2set)):
        tmp = []
        for t in range(len(Tset)):
            tmp.append(GetM(dat[j1j2,:,t,3],dat[j1j2,:,t,4],dat[j1j2,:,t,2]))
        datM0.append( np.array(tmp))
        
    datM0 = np.array(datM0).reshape((len(datM0),len(Tset),1))
    dath0 = np.concatenate((dath0,datM0),axis=2)
    
    #rawR.append(dath0)
    ## number of J1J2 sets
    for a in range(len(dath0)):
        f = open(os.path.join(Mdir,fname+'.%d.result'%(a)),'w')
        for t in range(len(dath0[a])):
            f.write("%11.11lf %11.11lf %11.11lf %11.11lf %11.11lf %11.11lf %11.11lf\n"%(dath0[a,t,0],dath0[a,t,1],dath0[a,t,2],dath0[a,t,3],dath0[a,t,4],dath0[a,t,5],dath0[a,t,6]))
        f.close()
    
    ## Save:
    


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
