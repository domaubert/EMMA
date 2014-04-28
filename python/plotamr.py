#!/usr/bin/env python
import sys, os
import numpy as np
import matplotlib.pylab as plt

from fonction_IO  import *
from fonction_amr import *
from plotpart import *



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

	args.field = ["x"]
	denoct2grid(filename, args,0)
	data = cube(filename + ".cube").getData()

	data =  np.sum( data,axis=0)
	data = np.log10(data)

	plt.imshow(data, interpolation='nearest')
	plt.colorbar()

	N,t, parts=ReadStars(args.files[0], args)
	plotpart(N,parts)


	plt.show()



def plotTR(args):
	rho = getField(args,args.files[0], "field.d")
	T   = getField(args,args.files[0], "temp")

	plt.loglog(rho/np.mean(rho),T, '.')


def plotStarTR(args):

	data = readFile("RT")

	data = [x for x in data.split('\n')]

	for i in range(len(data)):
		data[i] = [x for x in data[i].split('\t')] 

	n    = np.zeros(len(data))
	rho  = np.zeros(len(data))
	temp = np.zeros(len(data))
	z = np.zeros(len(data))
	ns=0

	for i in range(len(data)-1):
		if ( len(data[i])==6 and isfloat(data[i][1]) and isfloat(data[i][2]) and isfloat(data[i][3])  ) :
			if float(data[i][1]):
				ns +=1
				n[i] = float(data[i][1])
				rho[i] = float(data[i][2])/0.154330706937
				temp[i] = float(data[i][3])
				z[i] = float(data[i][5])
	

	lim1 = 3
	lim2 = 2

	mask = np.where(z>lim1)
	plt.loglog(rho[mask], temp[mask], 'r.')
	mask = np.where(np.logical_and(z<lim1, z>lim2) )
	plt.loglog(rho[mask], temp[mask], 'y.')
	mask = np.where(z<lim2)
	plt.loglog(rho[mask], temp[mask], 'g.')



def plotdiag(args):

	plotTR(args)
	plotStarTR(args)

	plt.xlim(1e-2, 1e4)
	plt.ylim(1e-2, 1e8)
	
	#plt.legend(getZ(args.files[0]))

	plt.show()




if __name__ == "__main__":	


	args = getargs()

#	oct2silo(args)
#	plothisto(args, field)
	plotslice(args)
#	plotdiag(args)


