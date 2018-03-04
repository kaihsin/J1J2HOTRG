import numpy as np
import torch
import os,sys
import math
from Util import *
import torch.nn.functional as F

"""
class Param:
    def __init__(self):
        self.J1   = None
        self.J2   = None
        self.Beta = None
        self.chi  = None
        self.nL   = None
"""


def Launch(J1,J2,Beta,h,chi,nL,savPath):

    #A = IsingLocal(Beta)

    #print(A)
    A = MakeLocalv2(J1,J2,h,Beta)
    # A.dim=[4 x 4 x 4 x 4] (l,u,r,d)
    #print(A)
    #exit(1)


    A = torch.from_numpy(A)
    #print(A)

    logNrm = []
    Ash = A.size()
    tr = A.contiguous().view(Ash[0]*Ash[1],-1).trace()
    A /= tr
    logNrm.append(np.log(tr))
    print(tr)


    for i in range(nL):
        L = 2**(i+1)

        #Update LR:
        A = Update(A,chi,0)
        # Update UD:
        A = Update(A,chi,1)

        Ash = A.size()
        tr = A.contiguous().view(Ash[0]*Ash[1],-1).trace()
        A /= tr
        logNrm.append(np.log(tr))

        lnZ = 0

        for j in range(len(logNrm)):
            lnZ += logNrm[len(logNrm)-j-1] * 4**j / L**2
            #print(lnZ)

        F = -lnZ/Beta

        # Write to file
        f = open(os.path.join(savPath,"%d.dat"%(L)),'a+')
        f.write("%11.11lf %11.11lf %11.11lf %11.11lf %11.11lf %11.11lf\n"%(J1,J2,1./Beta,h,lnZ,F))
        f.close()
        print("L=%d J1=%4.6lf J2=%4.6lf T=%4.6lf h=%4.6lf lnZ=%4.6lf F=%4.6lf"%(L,J1,J2,1./Beta,h,lnZ,F))


#===============================
if not os.path.exists("Data"):
    os.system("mkdir Data")

ID = "J1_10_J2_00_Ky8v2"

J1 = -1.
J2 = 0.0

chi = 8
nL = 16


Ti = 2.0
Tf = 3.0
NT = 64

#hi  = -0.01
#hf  = 0.01
#Nh = 16
hi = 0
hf = 0 
Nh = 1

#prepare Path:
savDir = "Data/%s"%(ID)
if os.path.exists(savDir):
    os.system("rm -r %s"%(savDir))
os.system("mkdir %s"%(savDir))


for t in range(NT):
    for n in range(Nh):
        Beta = 1./(Ti + t*(Tf-Ti)/NT)
        h    = hi + n*(hf-hi)/Nh
        Launch(J1,J2,Beta,h,chi,nL,savDir)
    #exit(1)
