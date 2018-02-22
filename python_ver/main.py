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

    #A = IsingLocal(Beta)
    
    #print(A)
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
        f.write("%11.11lf %11.11lf %11.11lf %11.11lf %11.11lf\n"%(J1,J2,Beta,lnZ,F))
        f.close()
        print("L=%d J1=%4.6lf J2=%4.6lf Beta=%4.6lf lnZ=%4.6lf F=%4.6lf"%(L,J1,J2,Beta,lnZ,F))


#===============================	

if not os.path.exists("Data"):
	os.system("mkdir Data")

ID = "J1_10_J2_00_Ky12"

J1 = -1.
J2 = 0.0

chi = 12
nL = 16

Ti = 1.0
Tf = 3.0
NT = 128

#prepare Path:
savDir = "Data/%s"%(ID)
if os.path.exists(savDir):
	os.system("rm -r %s"%(savDir))
os.system("mkdir %s"%(savDir))


for t in range(NT):
    Beta = 1./(Ti + t*(Tf-Ti)/NT)	
    Launch(J1,J2,Beta,chi,nL,savDir)
    #exit(1)


