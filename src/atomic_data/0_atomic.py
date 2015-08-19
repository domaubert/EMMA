#!/usr/bin/env python
import numpy as np
from scipy import integrate

#localisation of the spectrum file
filename='./0_Starburst_spectrum.dat'

#frequency bounds
SPACE_BOUND = [13.6, np.inf];
#SPACE_BOUND = [13.6, 24.587387, 54.417760, np.inf];

# temporal bounds
TIME_BOUND=[0,np.inf]
#TIME_BOUND=[0,10,100,np.inf]

class Constantes : 
    """
    some physical constants
    """
    def __init__(self):

        self.c      = 299792458        # velocity of light in m/s
        self.planck = 6.62617e-34      # J/s Planck constant
        self.M0     = 1.9891e30        # Solar masse
        self.eV     = 1.602176565e-19  # electron Volt in Joules

def lambda2E(l):
    """
    convert wavelenght[m] to energy [J]
    """
    cst=Constantes() 
    c=cst.c
    h=cst.planck
    return h*c/l

def eV2lambda(eV):
    """
    convert energy [eV] to wavelenght [m]
    """
    cst=Constantes()
    h = cst.planck
    c = cst.c
    Ej = eV*cst.eV
    return h*c/Ej

def cross_section(egy):
    """
    compute cross section
    """
    P=2.963
    x=egy/4.298e-1
  
    a=(x-1.)**2
    b=np.power(x,P/2.-5.5)
    c=np.power(1+np.sqrt(x/3.288e1),-P)
    
    return 5.475e4*a*b*c*1e-22*(egy >= 13.6) # m2



#get constant
c = Constantes()

#reading header
file = open(filename, 'r')
file.readline()
file.readline()
name=file.readline().split()[2:]
file.close()
age = [	float(n[:-3]) for n in name]

#loading data
data = np.loadtxt(filename,unpack=True, skiprows=3)

#converting in SI
X = data[0]*1e-10 #A->m
data[1:]=np.power(10,data[1:]) *1e-7/1e-10 #log(Erg/A) -> J/m/s


NGRP_SPACE=len(SPACE_BOUND)-1
NGRP_TIME=len(TIME_BOUND)-1
NGRP = NGRP_SPACE*NGRP_TIME

print "NGRP_SPACE\t%d"%NGRP_SPACE

print "SPACE_BOUND(eV)  ",
for time in range(NGRP_TIME):
    print "%f  "%(SPACE_BOUND[time]),
print ""

print "NGRP_TIME\t%d"%NGRP_TIME

print "TIME_BOUND(MYR)  ",
for time in range(NGRP_TIME):
    print "%f  "%(TIME_BOUND[time]),
print ""
print ""

print "hnu            alphae             alphai             factgrp"

for time in range(NGRP_TIME):
    Epho = np.zeros( (NGRP, len(name)) )
    s_N = np.zeros( (NGRP, len(name)) )
    s_E = np.zeros( (NGRP, len(name)) )
    nphot = np.zeros( (NGRP, len(name)) )

    b1= TIME_BOUND[time]
    b2= TIME_BOUND[time+1]
    print ""

    for num in range(1,len(name)+1):
        if( (age[num-1]>b1) & (age[num-1]<=b2)):
            for igrp in range(NGRP_SPACE):
                #selecting spectrum
                Y=data[num]

                #keep only group photons
                b1= eV2lambda(SPACE_BOUND[igrp  ])
                b2= eV2lambda(SPACE_BOUND[igrp+1])

                mask=np.where( (X<=b1) & (X>b2) )
                y=Y[mask]
                x=X[mask]

                #integrate total ionising energy
                Eion= integrate.trapz(y,x) #J/s

                #computing mean photon energy
                Nphot_per_sec = integrate.trapz(y/lambda2E(x),x)
                Epho[igrp, num-1]= Eion/Nphot_per_sec

                #computing nphot/sec/kg
                nphot[igrp,num-1]=Nphot_per_sec/(1e6*c.M0)

                #compute photoionisation cross section
                s_lambda = cross_section(lambda2E(x)/c.eV)

                #compute sigma N
                N_lambda = y/lambda2E(x)
                s_N[igrp,num-1] = integrate.trapz(s_lambda*N_lambda,x) /Nphot_per_sec

                #compute sigma E
                s_E[igrp,num-1] = integrate.trapz(s_lambda*y,x)/Eion


    for igrp in range(NGRP_SPACE):
        Nphot=integrate.trapz(nphot[igrp],age)
        E_tot = integrate.trapz(Epho[igrp]*nphot[igrp],age)/ Nphot
        s_N_tot = integrate.trapz(s_N[igrp]*nphot[igrp],age)/ Nphot
        s_E_tot = integrate.trapz(s_E[igrp]*nphot[igrp],age)/ Nphot	
        factgrp = np.sum(nphot[igrp])/np.sum(nphot)

        print "%s  %s  %s  %s"%(E_tot/c.eV,s_N_tot,s_E_tot, factgrp)

