#!/usr/bin/env python

import time, sys
import numpy as np
import matplotlib.pylab as plt
from fonction_part import *
from fonction_physique import *




def spectre(N,t,parts) :


	M = m2mo(parts.mass,1)

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

	file = args.folder[0] + "part." + str(args.snap[0]).zfill(5)

	N,t, parts=readPart(file, args)		

#	plotpart(args,N,t,parts)
	spectre(N,t,parts)	

	plt.show()





