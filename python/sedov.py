#!/usr/bin/env python
import sys, os
import numpy as np
import matplotlib.pylab as plt

from fonction_IO  import *
from fonction_amr import *
from plotpart import *


def getProfile(fileName,args):

	denoct2grid(fileName, args,0)
	data =  cube(args.files[0] + ".cube")

	c =  data.getData()
	Nx,Ny,Nz = data.getSize()

	x = np.arange(Nx/2)
	y = c[Nx/2:Nx,Ny/2,Nz/2]

	return Nx,Ny,Nz,x,y


def getProfileDiag(fileName,args):

	denoct2grid(fileName, args,0)
	data =  cube(args.files[0] + ".cube")

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


if __name__ == "__main__":	


	args = getargs()

	fileName = args.files[0] 
	fileName = fileName[:-10] +  "grid" +  fileName[-6:] 

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


	args.field = ["field.u"]
	Nx,Ny,Nz,x,y = getProfileDiag(fileName,args)
	plt.plot(x,y, label="diag")

	Nx,Ny,Nz,x,y = getProfile(fileName,args)
	plt.plot(x,y, label="x")


	plt.legend()
	plt.show()

	
