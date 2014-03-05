import sys
import os
import numpy as np
import matplotlib.pylab as plt

from fonction_part import *
from fonction_amr import *




def checkMtot(filename):
	a=array(filename)
	data = a.getData()
	N = a.getN()
	aexp = a.geta()
	print N

	lmax=6
	vcell = pow(2,-3* lmax)
	
	for i in range(N):
		data[i]	 *= vcell 

	print data.sum()




if __name__ == "__main__":	

	filename = "../utils/den.00060"

	checkMtot(filename)
