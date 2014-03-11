import sys
import os
import numpy as np
import matplotlib.pylab as plt

from fonction_amr import *


def plothisto(filename):
	
	a=array(filename)
#	data =  np.log10(a.getData())
	data =  a.getData()
	n, bins, patches = plt.hist(data, 100, normed=0, histtype='stepfilled', log=True)
	plt.show()


def plotslice(filename):

	a=cube(filename)

	data = np.mean( a.getData(),axis=2)
	print np.max(a.getData())
#	data =  a.getData()[32,:,:]

	plt.imshow(data, interpolation='nearest')
	plt.colorbar()
	plt.show()





if __name__ == "__main__":	

	filename = "../utils/den.00060"

	plothisto(filename)
#	plotslice(filename)
 
