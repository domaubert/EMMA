#!/usr/bin/env python
from fonction_part import *
from fonction_IO import *
from fonction_physique import *
import numpy as np
 


def readTdata(folder):


	fname = folder + "tdata.00000.p00000"
	print fname
	fs = open(fname)	
	n = (int)(fs.read(8))
	N = np.fromfile(fs, dtype=np.int32, count=n) 
	A = np.fromfile(fs, dtype=np.float64, count=n) 
	Mtot = np.fromfile(fs, dtype=np.float64, count=n) 
	fs.close()
	return n,N,A,Mtot


if __name__ == "__main__":


	args=getargs()
	
	L = [10,10,10,2,2,2]
	Fol=0

	for folder in args.folder:

		n,N,A,Mtot = readTdata(folder)

		npart0 = N[0]
	
		for i in range(n) :
			N[i] -= npart0

		Z = np.zeros(n, dtype = np.float64)
		T = np.zeros(n, dtype = np.float64)
		SFR = np.zeros(len(Mtot))


		for i in range(n) :
			Z[i] = a2z(A[i])
			T[i] = a2t(A[i]) 



		for i in range(n-1) :
			SFR[i] = m2mo(Mtot[i+1] - Mtot[i], A[i],L[Fol]) / (T[i+1] - T[i]) / pow(L[Fol],3)
		Fol +=1
		f= findfirststar(N)+1
		plt.semilogy(Z[f:-1],SFR[f:-1], label = folder)


		



	plt.xlabel(r'$Z$')
	plt.ylabel(r'$SFR (M_{o}.yr^{-1}.Mpc^{-3})$' )
	plt.legend()

	plt.show()











