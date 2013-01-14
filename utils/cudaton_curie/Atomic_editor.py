#!/usr/bin/env python

import sys 
import os
from math import *
import numpy as np
from scipy.io.numpyio import fwrite,fread
from scipy.integrate import *

def begy(ee,T,s_t):
	k=1.3806e-23		#J/K1
	c=2.99792458e8		#m/s
	h=2*pi*1.054572e-34	#J.s
	
	if (s_t==1):
		ee=ee*1.6022e-19
		fac=2*pi/h/h/h/c/c

		aa=exp(ee/k/T)
		res=fac*pow(ee,3)/(-1+aa)	# m-2 s-1
	else:
		nu=ee/h/c
		res= pow (nu, -T)
	return res

def so(egy):
	sig0=5.475e4
	E0=4.298e-1
	y0=0.
	y1=0.
	yw=0.
	ya=3.288e1
	P=2.963
	egyth=13.6

	x=egy/E0-y0
	y=sqrt(x*x+y1*y1)

	F=( pow(x-1.,2) + yw*yw ) * pow(y,0.5*P-5.5) * pow( 1+sqrt(y/ya),-P)
	
	res=sig0*F*1e-22*(egy >= egyth)	# m2

	return res

def getsrcparam(mini,maxi):
	avg_se=np.zeros(len(mini))
	avg_sn=np.zeros(len(mini))
	avg_ee=np.zeros(len(mini))
	avg_en=np.zeros(len(mini))
	fact=np.zeros(len(mini))
	E=np.zeros(len(mini))
	N=np.zeros(len(mini))
	ev=1.6022e-19	

	T=[50000,100000,1.8,5]
	Tstr=["S_50000", "S_100000", "S_1p8", "S_5p0"]

	atomic=open("Atomic.h","wb")
	n="#define NGRP " + str(len(mini)-1) + "\n\n"
	atomic.write(n)

	for j in range(len(T)):
		s_t =  j<2

		for i in range(len(mini)):
			int_A = quad(lambda x : so(x)*begy(x,T[j],s_t),mini[i],maxi[i])
			int_B = quad(lambda x: begy(x,T[j],s_t),mini[i],maxi[i])
			avg_se[i] = int_A[0]/int_B[0]
	
			int_A = quad(lambda x : so(x)*begy(x,T[j],s_t)/x,mini[i],maxi[i])
			int_B = quad(lambda x: begy(x,T[j],s_t)/x,mini[i],maxi[i])
			avg_sn[i] = int_A[0]/int_B[0]
	
			int_A = quad(lambda x : begy(x,T[j],s_t)*x,mini[i],maxi[i])
			int_B = quad(lambda x: begy(x,T[j],s_t),mini[i],maxi[i])
			avg_ee[i] = int_A[0]/int_B[0]
	
			int_A = quad(lambda x : begy(x,T[j],s_t),mini[i],maxi[i])
			int_B = quad(lambda x: begy(x,T[j],s_t)/x,mini[i],maxi[i])
			avg_en[i] = int_A[0]/int_B[0]	
			N[i]=int_B[0]
		
		s0="#ifdef " + str(Tstr[j]) + "\n"
		atomic.write(s0)
		atomic.write("#define SECTION_EFFICACE ")
		for i in range(0,len(mini)-1):
			s1 = "hnu[" + str(i) + "]=" + str(avg_en[i]) + "*" + str(ev) + ";" 
			atomic.write(s1)
			s2 = "alphae[" + str(i) + "]=" + str(avg_se[i]) + "*c;" 
			atomic.write(s2)
			s3 = "alphai[" + str(i) + "]=" + str(avg_sn[i]) + "*c;"
			atomic.write(s3)
	
		atomic.write("\n#define FACTGRP ")
		for i in range(0,len(fact)-1):
			s4 = "factgrp[" + str(i) + "]=" + str(N[i]/N[len(fact)-1]) + ";"		
			atomic.write(s4)

		atomic.write("\n#endif\n\n")
	atomic.close()

def main():

	test=1
	while test :

		Def= raw_input("Utiliser les valeurs par defaut ? o/n : \n")
	
		if Def == 'o' or Def =='O' :
	
			mini=[13.6, 24.6, 54.4, 13.6]
			maxi=[24.6, 54.4, 1e3 , 1e3 ]

			test=0
	
		elif Def == 'n' or Def =='N' :
		
			Nb = raw_input("Nombre de Groupes : \n")
			Nb=int(Nb)
	
			mini=np.zeros(Nb+1)
			maxi=np.zeros(Nb+1)
	
			mini[0]=13.6
	
			for i in range (1,Nb):
				s="Fin du groupe N "+str(i) +"\n"
				maxi[i-1]=raw_input(s)
				mini[i]=maxi[i-1]
			maxi[Nb-1]=1e3
			mini[Nb]=13.6
			maxi[Nb]=1e3

			test=0	
		else :
			print " Entree incorrecte"
	
	G=getsrcparam(mini,maxi)

	print "Ecriture de 'Atomic.h' OK"

if __name__ == '__main__':
	main()


