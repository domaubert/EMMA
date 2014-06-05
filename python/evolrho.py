#!/usr/bin/env python
from fonction_part import *
from fonction_amr import *
from fonction_IO import *
from fonction_physique import *
import numpy as np
import matplotlib.pylab as plt 


def getThresh0(A):
	H0 		= 67 *1000.0/1e6/3.085677e16;
	rhostar 	=   3.0 * pow(H0,2.0) * 0.3175 /(8.0*np.pi*6.67384e-11);

	rhocrittilde 	= 360 * 1.67262158e-27;

	rhocrit = np.zeros(len(A))
	i = 0
	for a in A :
		rhocrit[i] 	= pow(a,3)* rhocrittilde / rhostar; 
		i+=1

	return rhocrit



def getThresh1(A):
	rhocrit = np.zeros(len(A))
	i = 0
	for a in A :
		rhocrit[i] 	= 0.154330706937 * 5
		i+=1

	return rhocrit



if __name__ == "__main__":
	args=getargs()
		
	n = len(args.files)
	N = pow(2.,3*args.level)

	Rho = np.zeros((n,N), dtype = np.float64)
	A   = np.zeros(n, dtype = np.float64)
	Z   = np.zeros(n, dtype = np.float64)

#	genCube(args, "field.d")

	i=0
	for file in args.files:
		A[i] = getA(file)
		Z[i] = getZ(file)
		print Z[i]
		Rho[i] = array(file[:-7] + ".cube").getData()
		i+=1


	
	thresh0=getThresh0(A)
	thresh1=getThresh1(A)
	plt.semilogy(Z,thresh0, linewidth=4)
	plt.semilogy(Z,thresh1, linewidth=4)



	snap = n-8
	print "Z = ",Z[snap]
	
	idx = np.where( Rho[snap][:] > thresh0[snap] )[0]
	print idx.shape



	tmp = np.zeros(n, dtype = np.float64)
	for i in idx:
		for file in range(n):
			tmp[file] = Rho[file][i]
		plt.semilogy(Z,tmp)





	plt.xlim(0,15)
#	plt.ylim(1e-3,10)

	plt.xlabel(r'$Z$')
	plt.ylabel(r'$Rho max$' )
	plt.legend(loc=1)

	plt.show()











