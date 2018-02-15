import numpy as np


def Derivative(x,y,fitorder=1,grppt=2):
	outx,outy = [],[]
	
	for i in range(len(y)-grppt):
		c = np.polyfit(x[i:i+grppt],y[i:i+grppt],deg=fitorder)
		cdiv = np.polyder(c)
		ix= np.mean(x[i:i+grppt])
		outy.append(np.polyval(cdiv,ix))
		outx.append(ix)
	return np.array(outx),np.array(outy)




