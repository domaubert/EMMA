#!/usr/bin/env python
import sys, os
import numpy as np
import matplotlib.pylab as plt
import math as ma

from fonction_IO  import *
from fonction_amr import *
from plotpart import *


def getProfile(fileName,args):
	level = args.level
	field = args.field[0]


	file  = fileName.replace("grid","cube"+ str(level)) + "." + field

	data =  cube(file)
	Nx,Ny,Nz = data.getSize()	
	c = data.getData()
	
	x = np.arange(Nx/2)
	y = c[Nx/2:Nx,Ny/2,Nz/2]
	
	return Nx,Ny,Nz,x,y

def getProfileDiag(fileName,args):
	level = args.level
	field = args.field[0]

	file  = fileName.replace("grid","cube"+ str(level)) + "." + field

	data =  cube(file)

	c =  data.getData()
	Nx,Ny,Nz = data.getSize()

#	x = np.arange(Nx/2)
#	y = c[Nx/2:Nx,Ny/2,Nz/2]

	diag = int(np.sqrt(Nx*Nx/4 + Ny*Ny/4))

	x = np.arange(Nx/2) * np.sqrt(2)
	y = np.zeros(Nx/2)
	for i in range(Nx/2):
		y[i] = c[Nx/2+i,Ny/2+i,Nz/2]

	return Nx,Ny,Nz,x,y

def getProfileProper(fileName,args):

	level = args.level
	field = args.field[0]

	file  = fileName.replace("grid","cube"+ str(level)) + "." + field

	data =  cube(file)

	c =  data.getData()
	Nx,Ny,Nz = data.getSize()


	xc,yc,zc = Nx/2.,Ny/2.,Nz/2.


	diag = int(Nx/2 * np.sqrt(3)) +1

	rho = np.zeros(diag, dtype=np.float)
	n = np.zeros(diag, dtype=np.float)
	xplot = np.arange(diag)

	for z in range(0,Nz):	
		print z
		Z = float(z-zc)
		for y in range(0,Ny):
			Y = float(y-yc)
			for x in range(0,Nx):
				X = float(x-xc)

				r = int(np.sqrt( X*X + Y*Y + Z*Z ) )
				rho[r] += c[z][y][x]
				n[r] += 1


	y = np.divide(rho,n)

	return Nx,Ny,Nz,xplot[:Nx/2],y[:Nx/2]



def solution():
	
	f  = open('python/sedov.dat','r')
	data = f.read()
	data = data.split('\n')[6:]
	
	n = len(data) -1
	x = np.zeros(n)
	y = np.zeros(n)
	
	
	for i in range(n):
		ligne = data[i].split('    ')
		x[i] =  float(ligne[1])
		y[i] =  float(ligne[2])

	f.close()

	return x,y
	

def solution2():
		return np.loadtxt('python/sedov2.dat', usecols=(1,2), unpack=True)

	
if __name__ == "__main__":	

	args = getargs()


	fileName = args.files[0] 
	fileName = fileName[:-10] +  "grid" +  fileName[-6:] 


	args.field = ["field.d"]
#	denoct2grid(fileName, args,0)



	"""
	args.field = ["field.d"]
	Nx,Ny,Nz,x,y = getProfile(fileName,args)
	plt.plot(x,y, label="density")

	args.field = ["field.p"]
	Nx,Ny,Nz,x,y = getProfile(fileName,args)
	plt.plot(x,y, label="P")

	args.field = ["field.u"]
	Nx,Ny,Nz,x,y = getProfile(fileName,args)
	plt.plot(x,y, label="Vx")
	"""

	"""
	args.field = ["field.d"]
	Nx,Ny,Nz,x,y = getProfile(fileName,args)
	plt.plot(x,y, label="x")

	Nx,Ny,Nz,x,y = getProfileDiag(fileName,args)
	plt.plot(x,y, label="diag")
	
	Nx,Ny,Nz,x,y = getProfileProper(fileName,args)
	plt.plot(x,y, label="proper")
	"""
	
	
	#Nx,Ny,Nz,x,y = getProfileProper(fileName,args)
	Nx,Ny,Nz,x,y = getProfile(fileName,args)
	plt.plot(x,y,'.', label="numerical")
	
	
		
	x,y = solution2()
	x /= np.max(x)
	x *= pow(2,args.level)/2
	
	y /= y[len(y)-1]
	
	#print x,y

	plt.plot(x,y,label='analytical')
	
	plt.xlim(0,64)
#	plt.ylim(0.975,1.01)

	

	plt.legend()
	plt.show()

	

