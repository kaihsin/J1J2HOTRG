import numpy as np 
import matplotlib.pyplot as plt 
import os,sys
Dir = '../Data/test'


fs  = [x for x in os.listdir(Dir) if not os.path.isdir(os.path.join(Dir,x)) and not 'result' in x]

Ls = []
Es = []
Cvs = [] 
Datls = []
for fn in fs:
    Ls.append(int(fn.split('.')[0].split('L')[-1]))
    path = os.path.join(Dir,fn)

    # open lnZ ,F 
    f = open(path,'r')
    lines = f.readlines()
    dat = []
    for line in lines:
        dat.append( np.array(line.strip().split(' '),dtype=np.float))
    dat = np.array(dat)
    Datls.append(dat)

    # open E
    f = open(path+'.resultE','r')
    lines = f.readlines()
    dat = []
    for line in lines:
        dat.append( np.array(line.strip().split(' '),dtype=np.float))
    dat = np.array(dat)
    Es.append(dat)

    # open Cv
    f = open(path+'.resultCv','r')
    lines = f.readlines()
    dat = []
    for line in lines:
        dat.append( np.array(line.strip().split(' '),dtype=np.float))
    dat = np.array(dat)
    Cvs.append(dat)



plt.figure(1)
# F
plt.title('F')
for l in range(len(Ls)):
    plt.plot(Datls[l][:,2],Datls[l][:,4],'o',label='L%d'%(Ls[l]))

plt.xlabel("T")
plt.ylabel("F")
plt.legend()
plt.grid(1)


plt.figure(2)
# E
plt.title('E')
for l in range(len(Ls)):
    plt.plot(Es[l][:,2],Es[l][:,3],'o',label='L%d'%(Ls[l]))

plt.xlabel("T")
plt.ylabel("E")
plt.legend()
plt.grid(1)

plt.figure(3)
# E
plt.title('Cv')
for l in range(len(Ls)):
    plt.plot(Cvs[l][:,2],Cvs[l][:,3],'o',label='L%d'%(Ls[l]))

plt.xlabel("T")
plt.ylabel("Cv")
plt.legend()
plt.grid(1)





