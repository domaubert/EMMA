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
	bkg0 = np.zeros(n)
	bkg1 = np.zeros(n)
	bkg2 = np.zeros(n)
	for z in range(n):
		bkg0[z] = amp/(sig*np.sqrt(2*np.pi))* np.exp(-pow((z-zavg),2)/(2.*pow(sig,2)))+mz*z+pz
		bkg0[z] *= 1e15

		if z<3 : 
			bkg1[z] = 3.6
		else :
			bkg1[z] = 3.6*(4./(1+z))

		if z>=2 : 
			bkg2[z] = 3.6*pow(1.+z,-3.)
		else :
			bkg2[z] = 3.6*(1+z)/81.




	plt.plot(range(n),bkg0,label = "Haardt & Madau")
	plt.plot(range(n),bkg1,label = "Katz")
	plt.plot(range(n),bkg2,label = "Moi")
	
def observation():
	Z   = [0.2,2,3.8,4.9,5.9,6.8,7.9,10.4]

	SFR1 = [-2.2,-1.4,-1.44,-1.53,-1.84,-1.95,-2.2,-3.33]
	SFR1 = [ pow(10,x) for x in SFR1]
	plt.plot(Z,SFR1, 'k--')
#, label = "obs dust uncorrected"

	
	SFR2 = [-1.7,-1,-1.06,-1.19,-1.59,-1.72,-2.05,-3.18]
	SFR2 = [ pow(10,x) for x in SFR2]
	plt.plot(Z,SFR2, 'k--')
#,label = "obs dust corrected"

def SFRFromTdata(args) :

	Fol=0


	for folder in args.folder:

		n,N,A,Mtot = readTdata(folder)

		npart0 = N[0]
		for i in range(n) :
			N[i] -= npart0

		Z = np.zeros(n, dtype = np.float64)
		T = np.zeros(n, dtype = np.float64)
		SFR = np.zeros(len(Mtot))

		for i in range(n-1) :
			Z[i] = a2z(A[i+1]) 
			T[i] = a2t(A[i+1]) 
		
		for i in range(n-1) :
			SFR[i] = ( m2mo(Mtot[i+1],1) - m2mo(Mtot[i],1) ) / (T[i+1] - T[i]) 
		print SFR


		for i in range(n-1) :
			Z[i] = (a2z(A[i+1]) + a2z(A[i]))/2.
			T[i] = (a2t(A[i+1]) + a2t(A[i]))/2.


		
		Fol +=1

		f= findfirststar(N)
		plt.semilogy(Z[f:],SFR[f:], label = folder) 					# Star Formation Rate en fonction du redshift
#		plt.plot(T[f:-1]/1e9,Mtot[f:-1]/0.154330706937 , label = folder)		# Masse d'etoile en fonction du temps




#	uvbkg()
	observation()

	plt.xlim(0,15)
#	plt.ylim(1e-3,10)
	plt.xlabel(r'$Z$')
	plt.ylabel(r'$SFR (M_{o}.yr^{-1}.Mpc^{-3})$' )
	plt.legend(loc=1)

	plt.show()

def SFRFromSnap(stars, lab):


	if len(stars.mass):
		b = 16
		n0, bin_edges = np.histogram(stars.age, bins = b)
	
		z  =a2z(t2a( bin_edges))

		n= len(n0)
		sfr = np.zeros(n)
		for i in range(n):
			sfr[i] = m2mo( stars.mass[0] * (n0[i]) ,1 )  / float( bin_edges[i+1] -bin_edges[i])
	
		plt.semilogy(z[:-1],sfr, label = lab)



if __name__ == "__main__":


	args=getargs()
	for file in args.files:
		Ntot,t,stars = readStars(file, args)
		SFRFromSnap(stars,file )



#	uvbkg()
	observation()

#	plt.xlim(0,15)
#	plt.ylim(1e-3,10)
#	plt.xlabel(r'$Z$')
#	plt.ylabel(r'$SFR (M_{o}.yr^{-1}.Mpc^{-3})$' )
	plt.legend()

	plt.show()











