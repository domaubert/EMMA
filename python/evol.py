#!/usr/bin/env python
from fonction_part import *
from fonction_IO import *
from fonction_physique import *
import numpy as np
import matplotlib.pylab as plt 
import matplotlib as mpl
from decimal import Decimal


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

	SFR2 = [-1.7,-1,-1.06,-1.19,-1.59,-1.72,-2.05,-3.18]
	SFR2 = [ pow(10,x) for x in SFR2]
	plt.semilogy(Z,SFR2, 'k-.',label = "obs dust corrected")

	SFR1 = [-2.2,-1.4,-1.44,-1.53,-1.84,-1.95,-2.2,-3.33]
	SFR1 = [ pow(10,x) for x in SFR1]
	plt.plot(Z,SFR1, 'k--', label = "obs dust uncorrected")

	



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

def findFirstStar(stars):
	return np.argmax(stars.age)

def plotFirtsFctRes(stars,P):
	
	id  = findFirstStar(stars)
	age = stars.age[id] 
	res = stars.mass[id] 


	size= int(P.unit_l/3.08567758e22 * 0.67 )
	dx = pow(2, -P.lmax)
#	res = size * dx
	
	
	label = "L = " + str(P.lcoarse) + "  --> " + str(size) + " Mpc/h"


	plt.plot(age,res, 'o', label=label)



def SFRFromSnap(stars, P, cool):

	if len(stars.mass):
		b = 8
		n0, bin_edges = np.histogram(stars.age, bins = b)

		n= len(n0)
		sfr = np.zeros(n)
		t   = np.zeros(n)
		for i in range(n):
			t[i]   = bin_edges[i] + (bin_edges[i+1] - bin_edges[i])/2 
			sfr[i] = m2mo( stars.mass[0] * (n0[i]) ,1)  / float( bin_edges[i+1] -bin_edges[i])

		z = a2z(t2a( t ))





		if P.lcoarse == 6:
			color = 'b'
		if P.lcoarse == 7:
			color = 'g'
		if P.lcoarse == 8:
			color = 'r'

		case = 10
		
		if case == 0 : # Etude de l'influence de delta
#  ./python/evol.py  -fi data/bu2/odc5/part.00001 -fi data/bu2/odc20/part.00001 -fi data/bu2/odc35/part.00001 -fi data/bu2/odc50/part.00001 -fi data/bu2/odc5_7/part.00001 -fi data/bu2/odc20_7/part.00001  -fi data/bu2/odc35_7/part.00001 -fi data/bu2/odc50_7/part.00001 -fi data/bu2/odc50_8/part.00001


			label = "L = " + str(P.lcoarse) + "   Delta = " + str(int(P.overdensity_cond))

			if P.overdensity_cond == 50:
				type = '-'
			if P.overdensity_cond == 35:
				type = '--'
			if P.overdensity_cond == 20:
				type = '-.'
			if P.overdensity_cond == 5:
				type = ':'		

		if case == 1 : # Etude de l'influence du raffinement 

#  ./python/evol.py -fi data/bu2/odc50/part.00001 -fi data/bu3/6_1/part.00001 -fi data/bu3/6_2/part.00001 -fi data/bu2/odc50_7/part.00001 -fi data/bu3/7_1/part.00001 -fi data/bu3/7_2/part.00001 -fi data/bu2/odc50_8/part.00001 

			dl  = P.lmax -  P.lcoarse
			label = "L = " + str(P.lcoarse) + " + " + str(dl)

			if dl == 0:
				type = '-'
			if dl == 1:
				type = '--'
			if dl == 2:
				type = '-.'
			if dl == 3:
				type = ':'
		
		if case == 2 : #Etude de l'influence des echelle spatiale

#  ./python/evol.py -fi data/bu4/12.5_6/part.00001 -fi data/bu4/25_7/part.00001 -fi data/bu2/odc50_8/part.00001 -fi data/bu4/25_6/part.00001 -fi data/bu2/odc50_7/part.00001 -fi data/bu4/100_8/part.00001 -fi data/bu2/odc50/part.00001 -fi data/bu4/100_7/part.00001 -fi data/bu4/200_8/part.00001 -fi data/bu4/12.5_6/part.00001


			size = P.unit_l/3.08567758e22 * 0.67
			size=  np.around(size,2)

			label = "L = " + str(P.lcoarse) + "  --> " + str(size) + " Mpc/h"
			dx = pow(2, -P.lmax)
			res = size * dx

			print size, dx, P.lcoarse, res, 50 * pow(2, -6)


			if res == 50 * pow(2, -6) :
				type = '-'
			if res == 50 * pow(2, -7) :
				type = '--'
			if res == 50 * pow(2, -8) :
				type = '-.'



		if case == 3 : #Etude de l'influence de la condition de  densite physique 
#  ./python/evol.py -fi data/bu2/odc5_7/part.00001 -fi data/bu7/7_odc5_10/part.00001 -fi data/bu7/7_odc5_100/part.00001 -fi data/bu7/7_odc5_250/part.00001 -fi data/bu7/7_odc5_500/part.00001  -fi data/bu7/7_odc5_750/part.00001

	


			label = " Rho = " + str(-P.density_cond)
			print P.density_cond
			if P.density_cond == 0 :
				type = '-'
			if P.density_cond == -10 :
				type = '--'
			if P.density_cond == -100 :
				type = '-.'
			if P.density_cond == -250 :
				type = ':'
			if P.density_cond == -500 :
				color = 'b'
				type = '-.'
			if P.density_cond == -750 :
				color = 'b'
				type = ':'

		if case == 4 : #Etude de l'infulence de feedback

#./python/evol.py -fi data/bu2/odc50_7/part.00001 -fi data/bu5/7_0.5/part.00001 -fi data/bu5/7_1/part.00001  -fi data/bu5/7_10/part.00001 -fi data/bu5/7_0.5_sanscool/part.00001 -fi data/bu5/7_1_sanscool/part.00001 -fi data/bu5/7_10_sanscool/part.00001 


			label = " Cooling = " + str(int(cool)) + " Feedback = " + str(P.feedback_eff)
	

			if cool : 
				color = 'b'
			else :
				color = 'r'
			
			if P.feedback_eff == 0:
				type = '-'
				color = 'g'
			if P.feedback_eff == 0.5:
				type = '--'
			if P.feedback_eff == 1:
				type = '-.'
			if P.feedback_eff == 10:
				type = ':'
			
		if case == 5 : # Etude de l'influence de l'intensite des sources
#./python/evol.py -fi data/bu2/odc50_7/part.00001 -fi data/bu6/7_50_1e10/part.00001 -fi data/bu6/7_50_1e20/part.00001 -fi data/bu6/7_50_1e30/part.00001 -fi data/bu6/7_50_1e15/part.00001                    


			label = "Src int = %.0e" % P.srcint

			if P.srcint == 0:
				type = '-'
			if P.srcint == 1e10:
				type = '--'
			if P.srcint == 1e15:
				type = '-.'
			if P.srcint == 1e20:
				type = ':'
			if P.srcint == 1e30:
				type = ':'
			

		if case == 6 : # Etude du temps caracteristique de formation
# ./python/evol.py -fi data/bu2/odc5_7/part.00001 -fi data/bu8/7_ocd5_9e9/part.00001 -fi data/bu8/7_ocd5_15e9/part.00001
		
			label = "T_car = %.1e" % P.tcar

			if P.tcar == 3e9:
				type = '-'
			if P.tcar == 9e9:
				type = '--'
			if P.tcar == 15e9:
				type = '-.'
			if P.tcar == 21e9:
				type = ':'

		if case == 7 : # etude vitesse de la lumiere 

			label = "C = %.1e" % P.clight
			
			if P.clight == 1e-2:
				type = '-'
			if P.clight == 1e-3:
				type = '--'
			if P.clight == 1e-4:
				type = '-.'			

		if case == 10 :
			label ="d = %d, rho = %.0e, Src = %.0e, Tcar = %.0e " % (P.overdensity_cond, P.density_cond, P.srcint, P.tcar)
			type = '-'


		plt.semilogy(z,sfr, linestyle=type, color=color, label = label)

		plt.xlabel(r'$Z$')
		plt.ylabel(r'$SFR$' )


def getXion(args):
	slurm = getSlurm(args)

	xion = []
	e = []
	z=  []
	print slurm
	for line in open(slurm):
		if "XION" in line:
			line = line.split('\t')
			xion.append(float(line[1]))
			z.append(line[3])

	plt.plot(z,xion, label="1-X")


def getE(args):
	slurm = getSlurm(args)

	e = []
	z=  []
	print slurm
	for line in open(slurm):
		if "EDEN\t" in line:
			line = line.split('\t')
			print line
			e.append(float(line[1]))
			z.append(line[3])

	plt.semilogy(z,e, label="1-X")

if __name__ == "__main__":

	args=getargs()
#	getXion(args)
	getE(args)

	cool = 0
	i    = 0
	for file in args.files:

		if i > 3 :
			cool = 1

		Ntot,t,stars = readStars(file, args)
	#	SFRFromSnap(stars,Param(file[:-10]),  cool)
		i += 1
#		plotFirtsFctRes(stars,Param(file[:-10]) )


	

#	uvbkg()
#	observation()


#	plt.title(r"L = 50 Mpc.h^-1")


	plt.xlim(4,11)
	plt.ylim(1e78, 1e85)
	plt.legend()
	plt.show()











