import numpy as np
import torch  
import os,sys
from Util import *
import torch.nn.functional as F

J1 = 1.
J2 = 0.
Beta = 1./0.2475
chi = 2

A = MakeLocal(J1,J2,Beta)
# A.dim=[2 x 2 x 2 x 2] (l,u,r,d)

A = torch.from_numpy(A)
print(A)
Ash = A.size()
tr = A.contiguous().view(Ash[0]*Ash[1],-1).trace()
A /= tr



for i in range(1):
    #Update LR:
    A = Update(A,chi)
    #Ash = A.size()
    #tr = A.contiguous().view(Ash[0]*Ash[1],-1).trace()
    #print(tr)
    #A /=tr
    # Update UD:
    A = Update(A,chi)
    


    Ash = A.size()
    tr = A.contiguous().view(Ash[0]*Ash[1],-1).trace()
    A /= tr
    print(tr)





