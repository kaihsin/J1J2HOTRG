import numpy as np 
import torch
import scipy as sp
from scipy import linalg


def StableSvd(M):

    #M is a torch tensor.
    npM = M.numpy()
    Q,R,P = sp.linalg.qr(npM,pivoting=True)    
    P = np.diag(np.ones(len(P)))[P,:]
    R = torch.from_numpy(R.dot(P))
    U,S,V = torch.svd(R)
    U = torch.matmul(torch.from_numpy(Q),U)

    return U,S,V

def MakeW(Beta):
    W = np.zeros((2,2))

    W[0,0] =  np.sqrt(np.cosh(Beta))
    W[0,1] =  np.sqrt(np.sinh(Beta))
    W[1,0] =  np.sqrt(np.cosh(Beta))
    W[1,1] = -np.sqrt(np.sinh(Beta))
    return W



def IsingLocal(Beta):
    tmp = np.zeros((2,2,2,2))
    W = MakeW(Beta)
    #print(W)
    #exit(1) 
    for l in range(2):
        for u in range(2):
            for r in range(2):
                for d in range(2):	
                    tmp[l,u,r,d] = W[0,l] * W[0,u] * W[0,r] * W[0,d] + W[1,l] * W[1,u] * W[1,r] * W[1,d]					
    return tmp
    

def MakeLocal(J1,J2,Beta):
    tmp = np.zeros((2,2,2,2))
     
    for l in range(2):
        for u in range(2):
            for r in range(2):
                for d in range(2):
                    sigD1 = (2.*l -1.)*(2.*d-1.)    
                    sigD2 = (2.*l-1.)*(2.*u-1.)
                    tmp[l,u,r,d] = (1.+(2.*l-1.)*(2.*r-1.)*(2.*u-1.)*(2.*d-1.))/2\
                                        * np.exp(-J1*Beta*0.5*(2.*(l+u+r+d)-4.) - J2*Beta*(sigD1+sigD2));

    return tmp


def Update(A,chi,direct):
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
    #U1 , S1, V1 = StableSvd(A)

    A = A.view(-1,Ash[2]*Ash[3]*Ash[4]*Ash[5])
    U2 , S2, V2 = torch.svd(A)
    #U2 , S2, V2 = StableSvd(A)
   
    #print("S1")
    #print(S1)
    #print("S2")
    #print(S2)
    ## truncation:
    new_chi = Ash[4]*Ash[5]

    if new_chi > chi:
        
        if torch.sum(S1[chi:]**2) > torch.sum(S2[chi:]**2):
            R = U2[:,:chi]
            L = R.transpose(0,1)
        else:
            R = V1[:,:chi]
            L = R.transpose(0,1)
        
        #R = U2[:,:chi]
        #L = R.transpose(0,1)
        #print(L)


        new_chi  = chi
    
        A = torch.matmul(L,A)
        A = A.view(-1,Ash[4]*Ash[5])
        A = torch.matmul(A,R)

        ## permute back :
    
    A = A.view(new_chi,Ash[2],Ash[3],new_chi) #[D,L,R,U]
    #A = A.permute(1,3,2,0)
    if direct ==0:
        A = A.permute(3,2,0,1) #output a 90 deg rotated permutation
    else:
        A = A.permute(0,1,3,2)
    return A


