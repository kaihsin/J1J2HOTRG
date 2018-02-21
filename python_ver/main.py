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


def Launch(J1,J2,Beta,chi,nL,savPath):

    A = MakeLocal(J1,J2,Beta)
    # A.dim=[2 x 2 x 2 x 2] (l,u,r,d)

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
        A = Update(A,chi)
        # Update UD:
        A = Update(A,chi)
    
        Ash = A.size()
        tr = A.contiguous().view(Ash[0]*Ash[1],-1).trace()
        if tr > 1E2:    
            A /= tr
        else:
            tr = 1
        print(tr)
        logNrm.append(np.log(tr))
	
        #lnZ = np.log(A.contiguous().view(Ash[0]*Ash[1],-1).trace())
        #print(lnZ)
        #lnZ /= L**2
        lnZ = 0	
        #print(logNrm)

        for j in range(len(logNrm)):
            lnZ += logNrm[len(logNrm)-j-1] * 4**j / L**2
            #print(lnZ)

        F = -lnZ/Beta		

        # Write to file
        f = open(os.path.join(savPath,"%d.dat"%(L)),'a+')
        f.write("%11.11lf %11.11lf %11.11lf %11.11lf %11.11lf\n"%(J1,J2,Beta,lnZ,F))
        f.close()
        print("L=%d J1=%4.6lf J2=%4.6lf Beta=%4.6lf lnZ=%4.6lf F=%4.6lf"%(L,J1,J2,Beta,lnZ,F))


#===============================	

if not os.path.exists("Data"):
	os.system("mkdir Data")

ID = "J1-10-J2-02-Ky10"

J1 = 1.
J2 = 0.2

chi = 10
nL = 12

Ti = 0.1
Tf = 1.0
NT = 32

#prepare Path:
savDir = "Data/%s"%(ID)
if os.path.exists(savDir):
	os.system("rm -r %s"%(savDir))
os.system("mkdir %s"%(savDir))


for t in range(NT):
	Beta = 1./(Ti + t*(Tf-Ti)/NT)	
	Launch(J1,J2,Beta,chi,nL,savDir)



