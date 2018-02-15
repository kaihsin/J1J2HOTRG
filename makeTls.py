import os,sys

pth = "Tls"

NJ1 = 1
J1i = 1.0
J1f = 1.0

NJ2 = 1
J2i = 0.0
J2f = 0.0

NJT = 80
Ti  = 0.2
Tf  = 4.0



f = open(pth,"w")
for i in range(NJ1):
    for j in range(NJ2):
        for k in range(NJT):
            f.write("%11.11lf %11.11lf %11.11lf\n"%(\
                    J1i + i*(J1f - J1i)/NJ1,\
                    J2i + j*(J2f - J2i)/NJ2,\
                    Ti  + k*(Tf  - Ti )/NJT))

f.close()





