#!/usr/bin/env python

import numpy as np
import matplotlib.pylab as plt
from fonction_part import *
from fonction_physique import *


def plotpart(N,parts):
	N=len(parts)
	x=np.zeros(N)
	y=np.zeros(N)

	for i in range(N) :
		x[i]=parts[i].x
		y[i]=parts[i].y

	Lmax=6
	N=pow(2,Lmax)
	H, xedges, yedges = np.histogram2d(x, y, bins=(N,N))

	plt.imshow(H, interpolation='nearest')
	plt.colorbar()


def spectre(N,t,parts) :

	M=[]
	for i in range (0,N):
		M.append(m2mo(parts[i].mass, t, 10))
		#M.append(parts[i].mass)

#	plt.hist(M, 100,  log=True, label = r'$Z = $' + str(a2z(t)).zfill(4) )	
	plt.hist(np.log10(M), 100,  log=True, label = r'$Z = $' + str(a2z(t)).zfill(4) )

	plt.legend()
#	plt.xlim(2,7)
#	plt.ylim(1e-1,1e5)
	plt.xlabel(r'$LOG 10 Mass [M_0]$')
	plt.ylabel(r'$PDF$' )
#	plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))



if __name__ == "__main__":
	
	args= getargs()

	N,t, parts=ReadStars(args.files[0], args)		

#	plotpart(N,parts)
	spectre(N,t,parts)	

	plt.show()





