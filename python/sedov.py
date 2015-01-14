#!/usr/bin/env python
import os, sys
import numpy as np
from pyemma  import *
import matplotlib.pylab as plt
from math import sqrt
import threading



def solution(level):
	"""
	eblast need to be set manually in sedov3.f 
	"""
	
	param = IO.Params().get()
	inyrs = float(param["unit_t"])/3.1536e13;
	
	a = amr.geta(amr.cubename(fileName,level,field))
	a -= 2e1
	a /= inyrs
	
	commande = "gfortran utils/sedov/sedov3.f -o utils/sedov/BW" 
	#os.system(commande)

	commande = "./utils/sedov/BW %s"%a
	os.system(commande)

	x,y = np.loadtxt('utils/sedov/sedov.dat', usecols=(1,2), unpack=True)
	x /= np.max(x)
	x *= pow(2,level)
	y /= y[len(y)-1]	
	plt.plot(x,y,label='analytical')
	
def plot_xft(args):
	data = readAllProfil(folder = args.folder)	
	t = amr.gettmap(folder = args.folder)
	unit = Cell2Meter(args)
	y = findFrontPosition(data)
	y = np.multiply(y,unit)

	plt.plot(t[1:],y[1:])
	plt.xlim(0,90)
	plt.ylim(0,600)
	plt.ylabel(r'position (pc)')
	plt.xlabel(r'time (Myr)')
		
def plot_vft(args):
	data = readAllProfil(folder = args.folder)	
	t = amr.gettmap(folder = args.folder)
	unit = Cell2Meter(args)

	y = findFrontPosition(data)
	y = np.multiply(y,unit)

	y2 = frontSpeed(y,t)
	
	
	plt.figure()
	plt.ylabel(r'speed (pc.Myr$^{-1})$')
	plt.xlabel(r'time (Myr)')
	plt.plot(t[1:],y2)
	
def plot_ampft(args):
	data = readAllProfil(folder = args.folder)	
	t = amr.gettmap(folder = args.folder)
	unit = Cell2Meter(args)

	y = findFrontAmp(data)

	plt.figure()
	plt.ylabel(r'amplitude')
	plt.xlabel(r'time (Myr)')
	plt.plot(t[1:],y[1:])
	plt.show()
	
if __name__ == "__main__":	

	args = IO.getargs()

#	Nx,Ny,Nz,x,y = getProfileProper(amr.cubename(args.files[0],args.level,args.field))
#	Nx,Ny,Nz,x,y = getProfile(fileName,args.level,args.field[0])
#	plt.plot(x,y,'.', label="numerical")
#	solution(args.level)

#	folders = ["data/", "src_sn/", "src/"]
	
#	for folder in folders:
#		args.folder=folder

#	getAllProfil(args, folder = args.folder)
#	plotAllProfil(args)

	plot_xft(args)
	#plot_vft(args)
	#plot_ampft(args)

	plt.show()
