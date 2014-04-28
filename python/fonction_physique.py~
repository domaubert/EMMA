import sys
import os 
import numpy as np

class Constantes : 
	def __init__(self):

		self.Parsec = 3.08567758e16
		self.G      = 6.67384e-11            	# newton constant
		self.MO     = 1.9891e30              	# Solar masse
		self.c      = 299792.458              	# velocity of light in km/sec

		self.Tyr = 977.8 

		self.H0 = 67.0			     	# Hubble constant 
		self.h  = self.H0/100.

		self.WM = 0.3175                 	# Omega(matter)
		self.WV = 0.6825               		# Omega(vacuum) or lambda
		self.WR = 4.165E-5/(self.h*self.h)    	# Omega(radiation)
		self.WK = 1-self.WM-self.WR-self.WV	# Omega curvaturve = 1-Omega(total)
	
def a2z(a) :
	return 1.0/a -1.0


def a2t(a) :
	c = Constantes()

	az = 1.0/(1.0+1.0*a2z(a))
	age = 0.
	n=1000         # number of points in integrals
	for i in range(n):
		a = az*(i+0.5)/n
		adot = np.sqrt(c.WK+(c.WM / a)+(c.WR/ (a*a) )+ (c.WV*a*a) )
		age = age + 1./adot
	zage = az*age/n
	zage_Gyr = (c.Tyr/c.H0)*zage

	return zage_Gyr*1e9



def m2mo(m, L) :
	c = Constantes()
	H0 = c.H0 * 1000.0/1e6/c.Parsec

	rho = 3.0*pow(H0,2.0) * c.WM / (8.0 * np.pi * c.G)
	
	L *= 1e6 * c.Parsec #Mpc en m
	V   = pow(L,3.0)

	Mtot = rho * V

	return Mtot * m / c.MO


def findfirststar(N) :
	first = 0
	while ( N[first] == 0 ) :
		first +=1
		if (first==len(N)) :
			break
	return first


