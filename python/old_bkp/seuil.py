#!/usr/bin/env python
from fonction_part import *
from fonction_IO import *
from fonction_physique import *
import numpy as np
import matplotlib.pylab as plt 
from matplotlib.ticker import NullFormatter

if __name__ == "__main__":
	PARSEC = 3.085677e16
	H0 = 67*1000.0/1e6/PARSEC;
	Om = 0.3175
	PROTON_MASS = 1.67262158e-27
	NEWTON_G  = 6.67384e-11
	aexp = 1 

	a3rhostar 	= pow(aexp,3.0) *  3.0 * pow(H0,2.0) * Om /(8.0*np.pi*NEWTON_G)

	
	rhocrittilde = 36 * PROTON_MASS
	rhocrit	= rhocrittilde / a3rhostar

	print rhocrit
