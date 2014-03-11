import sys
import os 
import numpy as np


def a2z(a) :
	return 1.0/a -1.0


def a2t(a) :

	G  = 6.67384e-11
	H0 = 67     # Hubble constant
	WM = 0.3175 # Omega(matter)
	WV = 0.6825 # Omega(vacuum) or lambda
	WR = 0.     # Omega(radiation)
	WK = 0.     # Omega curvaturve = 1-Omega(total)
	Tyr = 977.8 

	az = 1.0/(1+1.0*a2z(a))
	age = 0.
	n=1000         # number of points in integrals

	for i in range(n):
		a = az*(i+0.5)/n
		adot = np.sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
		age = age + 1./adot
	zage = az*age/n
	zage_Gyr = (Tyr/H0)*zage
	return zage_Gyr*1e9



def m2mo(m,a, L) :
	Parsec = 3.08567758e16
	G  = 6.67384e-11 #newton constant
	H0 = 67 * 1000/1e6/Parsec    # Hubble constant --> en si
	WM = 0.3175 # Omega(matter)
	WV = 0.6825 # Omega(vacuum) or lambda
	MO = 1.9891e30 #Solar masse

	rho = 3*pow(H0,2) *WM / (8 * np.pi * G)
	
	L *= 1e6 * Parsec #Mpc en m
	V   = pow(L,3)

	Mtot = rho * V

	return Mtot * m / MO

def findfirststar(N) :
	first = 0
	while ( N[first] == 0 ) :
		first +=1
		if (first==len(N)) :
			break
	return first


