import numpy as np
import IO

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

def z2a(z):
	return 1.0/(1.0+z)

def a2t(a) :
	"""
		convert expansion factor to time
	"""
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

def t2a(t):
	"""
		convert time to expansion factor
	"""
	n = 10000

	A = np.arange(n+1) / float(n)
	T = a2t(A)

	return np.interp(t,T,A)

def m2mo(m,folder) :
	"""
		convert code mass in kg
	"""
	c = Constantes()

	param = IO.Params(folder = folder).get()

	unit_m = float(param["unit_mass"])
	unit_m /= c.MO
	return m*unit_m 

def Cell2Meter(args):
	"""
		return the size of a cell in parsec
	"""
	param = IO.Params(folder = args.folder).get()
	L = float(param["unit_l"])/3.085677e16
	dx = pow(2,- args.level)*L
	return dx
