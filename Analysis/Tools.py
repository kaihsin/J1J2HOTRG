import numpy as np


def SmoothDiv(x,y,fitorder=1,grppt=2):
	outx,outy = [],[]
	
	for i in range(len(y)-grppt):
		c = np.polyfit(x[i:i+grppt],y[i:i+grppt],deg=fitorder)
		cdiv = np.polyder(c)
		ix= np.mean(x[i:i+grppt])
		outy.append(np.polyval(cdiv,ix))
		outx.append(ix)
	return np.array(outx),np.array(outy)


def GetE_smooth(Ts,lnZs):
    newT , dlnZs= SmoothDiv(Ts,lnZs,2,7) 
    return newT , newT**2 * dlnZs


def GetCv_smooth(Ts,lnZs):
    
    newTs,Es = GetE(Ts,lnZs)
    newT ,dEs = SmoothDiv(newTs, Es,1,3)
 
    return  newT, dEs


def SimpleDiv(x,y):
    newx = (x[1:] + x[0:-1])/2

    dydx = ( y[1:] - y[0:-1] ) / (x[1:] - x[0:-1]) 
    
    return newx,dydx


def GetE_simple(Ts,lnZs):
    newT, dlnZdT = SimpleDiv(Ts,lnZs) 
    return newT , newT**2 * dlnZdT


def GetCv_simple(Ts,lnZs):
    
    newTs,Es = GetE(Ts,lnZs)
    newT, dEdT = Div(newTs, Es)
 
    return  newT, dEdT



def GetM(hs,lnZs,Ts):
    print(hs)
    print(lnZs)
    print(Ts) 
    c  = np.polyfit(hs,lnZs*Ts,3)
    dc = np.polyder(c)
    
    return np.polyval(dc,0)





