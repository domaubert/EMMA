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



def uvbkg():

	amp=1.2e-16
	sig=1.
	zavg=2
	mz=1e-18
	pz=1.2e-17

	n=10
	bkg = np.zeros(n)
	for z in range(n):
		bkg[z] = amp/(sig*np.sqrt(2*np.pi))* np.exp(-pow((z-zavg),2)/(2.*pow(sig,2)))+mz*z+pz
		bkg[z] *= 1e15

	plt.plot(range(n),bkg,label = "UVbkg")


if __name__ == "__main__":


	args=getargs()
	
	Fol=0

#	h=62.0/100.0
#	L = [x * h for x in L] 

	for folder in args.folder:

		n,N,A,Mtot = readTdata(folder)
	#	print Mtot

	
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
			SFR[i] = ( m2mo(Mtot[i+1],1) - m2mo(Mtot[i],1) ) / (T[i+1] - T[i]) 
		print SFR


		Fol +=1

		f= findfirststar(N)
		plt.semilogy(Z[f:],SFR[f:], label = folder) 					# Star Formation Rate en fonction du redshift
#		plt.plot(T[f:-1]/1e9,Mtot[f:-1]/0.154330706937 , label = folder)		# Masse d'etoile en fonction du temps


	uvbkg()


	plt.xlim(0,15)
#	plt.ylim(1e-3,10)
	plt.xlabel(r'$Z$')
	plt.ylabel(r'$SFR (M_{o}.yr^{-1}.Mpc^{-3})$' )
	plt.legend(loc=1)

	plt.show()











