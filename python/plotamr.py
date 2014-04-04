#!/usr/bin/env python
import sys, os
import numpy as np
import matplotlib.pylab as plt

from fonction_IO  import *
from fonction_amr import *


def plothisto(args,field):
	
	filename = args.files[0][:-7]
	filename = filename[:-10] +  "grid" +  filename[-6:] 
	denoct2grid(filename, args, field,0)

	a=array("cube")
#	data =  np.log10(a.getData())
	data =  a.getData()
	n, bins, patches = plt.hist(data, 100, normed=0, histtype='stepfilled', log=True)
	plt.show()


def plotslice(args):

	filename = args.files[0][:-7]
	filename = filename[:-10] +  "grid" +  filename[-6:] 

	denoct2grid(filename, args,0)

	a=cube("cube")

	data = np.sum( a.getData(),axis=2)
	print np.max(a.getData())

#	data =  a.getData()[32,:,:]

	plt.imshow(np.log10(data), interpolation='nearest')
	plt.colorbar()
	plt.show()


def plotdiag(args):
	filename = args.files[0][:-7]
	filename = filename[:-10] +  "grid" +  filename[-6:] 

	denoct2grid(filename, args, "field.d",0)
	rho = array("cube").getData()

	denoct2grid(filename, args, "temp",0)
	T =array("cube").getData()

	print "Rho max", rho.max()
	print "T max", T.max()
	
	r = rho/np.mean(rho)
	print "Rho moyen", np.mean(rho)

	tmpR = []
	tmpT = []

	for i in range(len(r))	:
		if r[i]>2500 : 
			tmpR.append(r[i])			
			tmpT.append(T[i])			

	print len(tmpR)

	plt.loglog(r,T, '.')
	#plt.loglog(tmpR,tmpT, 'r.')
	plt.show()

def oct2silo(args):

	filename = args.files[0][:-7]
	filename = filename[:-10] +  "grid" +  filename[-6:] 
	denoct2grid(filename, args, args.field, 1)


if __name__ == "__main__":	


	args = getargs()
	field = "field.d"

#	oct2silo(args)
#	plothisto(args, field)
	plotslice(args)
#	plotdiag(args)


