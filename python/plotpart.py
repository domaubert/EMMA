#!/usr/bin/env python

import numpy as np
import matplotlib.pylab as plt
from fonction_part import *
from fonction_physique import *


def plotpart(N,parts):
	N=len(parts)
	x=np.zeros(N)
	y=np.zeros(N)

	Lmax=7
	S=pow(2.0,Lmax)

	for i in range(N) :
		x[i]=parts[i].x * S
		y[i]=parts[i].y * S

	

	'''
	H, xedges, yedges = np.histogram2d(x, y, bins=(N,N))
	plt.imshow(np.log10(H), interpolation='nearest')
	plt.colorbar()
	'''

	plt.plot(x,y,'.')

def spectre(N,t,parts) :

	M=[]
	for i in range (0,N):
		M.append(m2mo(parts[i].mass, t, 10))
		#M.append(parts[i].mass)

	lab =  r'$Z = $' + str(a2z(t)).zfill(4) 
#	plt.hist(M, 100,  log=True, label = lab )	
	plt.hist(np.log10(M), 100,  log=True, label = lab )

	plt.legend()
#	plt.xlim(2,7)
#	plt.ylim(1e-1,1e5)
	plt.xlabel(r'$LOG 10 Mass [M_0]$')
	plt.ylabel(r'$PDF$' )
#	plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))



if __name__ == "__main__":
	
	args= getargs()

	N,t, parts=ReadStars(args.files[0], args)		

	plotpart(N,parts)
#	spectre(N,t,parts)	

	plt.show()





