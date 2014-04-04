#!/usr/bin/env python
from fonction_part import *
from fonction_IO import *
from fonction_physique import *
import numpy as np
import matplotlib.pylab as plt 


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
	
	L = [10,20,20,2.5,5,40,2]
	Fol=0

	h=62.0/100.0

	L = [x * h for x in L] 


	for folder in args.folder:

		n,N,A,Mtot = readTdata(folder)
		print Mtot
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
			SFR[i] = ( m2mo(Mtot[i+1],L[Fol]) - m2mo(Mtot[i],L[Fol]) ) / (T[i+1] - T[i]) / pow(L[Fol],3.0)

		Fol +=1
		f= findfirststar(N)
		print SFR
		plt.semilogy(Z[f:-1],SFR[f:-1], label = folder) 				# Star Formation Rate en fonction du redshift
#		plt.plot(T[f:-1]/1e9,Mtot[f:-1]/0.154330706937 , label = folder)		# Masse d'etoile en fonction du temps


		


#	plt.xlim(0,25)
#	plt.ylim(1e-3,10)
	plt.xlabel(r'$Z$')
	plt.ylabel(r'$SFR (M_{o}.yr^{-1}.Mpc^{-3})$' )
	plt.legend(loc=1)

	plt.show()











