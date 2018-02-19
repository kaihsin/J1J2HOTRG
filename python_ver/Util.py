import numpy as np 
import torch
def MakeLocal(J1,J2,Beta):
    tmp = np.zeros((2,2,2,2))
     
    for l in range(2):
        for u in range(2):
            for r in range(2):
                for d in range(2):
                    sigD1 = (2.*l -1.)*(2.*r-1.)    
                    sigD2 = sigD1*(2.*l-1.)*(2.*u-1.)
                    tmp[l,u,r,d] = (1.+(2.*l-1.)*(2.*r-1.)*(2.*u-1.)*(2.*d-1.))/2\
                                        * np.exp(J1*Beta*0.5*(2.*(l+u+r+d)-4.) + J2*Beta*(sigD1+sigD2));

    return tmp


def Update(A,chi):
    ### Update LR
    # copy 
    B = A.clone()

    ## reshape A
    A = A.permute(0,1,3,2)
    Ash = A.size()
    A = A.contiguous().view(-1,Ash[3])

    ## reshape B
    Bsh = B.size()
    B = B.view(Bsh[0],-1)

    ## contract
    A = torch.matmul(A,B).view(Ash[0],Ash[1],Ash[2],Bsh[1],Bsh[2],Bsh[3])
    A = A.permute(2,5,0,4,1,3)
    Ash = A.size()

    A = A.contiguous().view(Ash[0]*Ash[1]*Ash[2]*Ash[3],-1)
    U1 , S1, V1 = torch.svd(A)
    A = A.view(-1,Ash[2]*Ash[3]*Ash[4]*Ash[5])
    U2 , S2, V2 = torch.svd(A)
    print("S1")
    print(S1)
    print("S2")
    print(S2)
    ## truncation:
    new_chi = Ash[4]*Ash[5]

    if new_chi > chi:
        V1 = V1[:,:chi]
        U2 = U2[:,:chi]
        new_chi  = chi
    
    U2.transpose_(0,1)

    A = torch.matmul(U2,A)
    A = A.view(-1,Ash[4]*Ash[5])
    A = torch.matmul(A,V1)

    ## permute back :
    A = A.view(new_chi,Ash[2],Ash[3],new_chi) #[D,L,R,U]
    #A = A.permute(1,3,2,0)
    A = A.permute(3,2,0,1) #output a 90 deg rotated permutation
    return A


