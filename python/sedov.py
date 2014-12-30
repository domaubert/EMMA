#!/usr/bin/env python
import os, sys
import numpy as np
from math import sqrt
import matplotlib.pylab as plt
import threading


from pyemma  import *



def getProfile(fileName,level,field):
	file = amr.cubename(fileName,level,field)
	data =  amr.cube(file)
	Nx,Ny,Nz = data.getSize()	
	c = data.getData()
	
	x = np.arange(Nx/2)
	y = c[Nx/2:Nx,Ny/2,Nz/2]
	
	return Nx,Ny,Nz,x,y

def getProfileDiag(fileName,args):
	level = args.level
	field = args.field[0]

	file  = fileName.replace("grid","cube"+ str(level)) + "." + field

	data =  amr.cube(file)

	c =  data.getData()
	Nx,Ny,Nz = data.getSize()

#	x = np.arange(Nx/2)
#	y = c[Nx/2:Nx,Ny/2,Nz/2]

	center = 0
	
	if center:
		diag = int(np.sqrt(Nx*Nx/4 + Ny*Ny/4))

		x = np.arange(Nx/2) * np.sqrt(2)
		y = np.zeros(Nx/2)
		for i in range(Nx/2):
			y[i] = c[Nx/2+i,Ny/2+i,Nz/2]
	else:
		diag = int(np.sqrt(Nx*Nx + Ny*Ny))

		x = np.arange(Nx) * np.sqrt(2)
		y = np.zeros(Nx)
		for i in range(Nx):
			y[i] = c[i,i,i]
		

	return Nx,Ny,Nz,x,y

class kernel(threading.Thread):
    def __init__(self, args):
        threading.Thread.__init__(self)
        self.args = args  
    def run(self):
		print "starting thread", self.threadID

   		c,rho,n,Nx,Ny,Nz,diag, Zmin,Zmax,xc,yc,zc, threadLock = self.args
   		
		for z in range(Zmin,Zmax):
			Z = float(z-zc)
			for y in range(0,Ny):
				Y = float(y-yc)
				for x in range(0,Nx):
					X = float(x-xc)

					r = int(np.sqrt( X*X + Y*Y + Z*Z ) )
					rho[r] += c[z][y][x]
					n[r] += 1
		threadLock.acquire()		
		threadLock.release()

def getProfileProperPara(fileName):

	data = amr.cube(fileName)

	Nx,Ny,Nz = data.getSize()
	xc,yc,zc = Nx/2.,Ny/2.,Nz/2.

	diag = int(Nx/2. * sqrt(3.0) ) +1

	nproc = 4

	deltaZ = Nz/nproc

	c =  data.getData()
	rho = np.zeros(diag, dtype=np.float)
	n = np.zeros(diag, dtype=np.float)
	
	thread=[]
	threadLock = threading.Lock()

	for i in range(nproc):			
		Zmin = deltaZ * (i  )
		Zmax = deltaZ * (i+1)
		
		thread.append(kernel((c,rho,n,Nx,Ny,Nz,diag,Zmin,Zmax,xc,yc,zc,threadLock)))
		
	for i in range(nproc):	
		thread[i].start()
		
	for i in range(nproc):	
		thread[i].join()
	print "all threads done"
	
	yplot = np.divide(rho,n)
	xplot = np.arange(diag)

	print np.sum(n), pow(256,3)

	return Nx,Ny,Nz,xplot,yplot

def getProfileProper(fileName, xc=0,yc=0,zc=0):

	data = amr.cube(fileName)

	Nx,Ny,Nz = data.getSize()
#	xc,yc,zc = Nx/2.,Ny/2.,Nz/2.
	
	diag = int((Nx-xc) * sqrt(3.0) ) +1

	c =  data.getData()
	rho = np.zeros(diag, dtype=np.float)
	n = np.zeros(diag, dtype=np.float)
	  
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
	yplot = np.divide(rho,n)
	xplot = np.arange(diag)

	return Nx,Ny,Nz,xplot,yplot
	
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

def comparField(args):
	args.field = ["field.d"]
	Nx,Ny,Nz,x,y = getProfile(fileName,args)
	plt.plot(x,y, label="density")

	args.field = ["field.p"]
	Nx,Ny,Nz,x,y = getProfile(fileName,args)
	plt.plot(x,y, label="P")

	args.field = ["field.u"]
	Nx,Ny,Nz,x,y = getProfile(fileName,args)
	plt.plot(x,y, label="Vx")
	
def comparMethod(args):
	args.field = ["field.d"]
	Nx,Ny,Nz,x,y = getProfile(fileName,args)
	plt.plot(x,y, label="x")
	
	Nx,Ny,Nz,x,y = getProfileDiag(fileName,args)
	plt.plot(x,y, label="diag")
	
	Nx,Ny,Nz,x,y = getProfileProper(fileName,args)
	plt.plot(x,y, label="proper")

def getAllProfil(args, folder = "data/", force=0):
	files = np.sort(os.listdir(folder))
	for file in files:
		if "grid." in file and ".p00000" in file: 
			if file.replace("grid","profile")[:-7] in files:
				print file.replace("grid","profile") + " allready exist"
				continue
			file= "%s%s"%(folder,file[:-7])

			convert.oct2grid(file,args.level,force=force)
			Nx,Ny,Nz,x,y = getProfileProper(amr.cubename(file,args.level,args.field[0]))				
			np.savetxt(file.replace("grid","profile"),(x,y))			

def gettmap(folder = "data/"):
	tmap = []
	for file in np.sort(os.listdir(folder)):
		if "grid." in file and ".p00000" in file: 
			file= "%s%s"%(folder,file)
			t = amr.geta(file)
			tmap.append(t)
	return np.array(tmap)

def readAllProfil(folder="data/"):
	data = []
	for file in np.sort(os.listdir(folder)):
		if "profile." in file : 
			data.append( np.loadtxt("%s%s"%(folder,file)) )
	return data

def plotAllProfil(args):
	data = readAllProfil(folder = args.folder)	

	for set in data:
		plt.plot(set[0],set[1])
	plt.show()

def findFrontPosition(data):
	position = np.zeros(len(data))
	i=0
	for set in data:
		grad = np.gradient(set[1])

		xinterp = np.linspace(0,len(data),4196)
		yinterp = np.interp(xinterp, set[0], grad)

		signchange = (np.diff(np.sign(yinterp))!=0)*1
		wheresignchange = np.where(signchange==1)
		
		if np.size(wheresignchange):			
			ind = np.argmax(set[1][np.int_(xinterp[wheresignchange[0]])]) #peut etre revoir l'arrondis ??
			pos = wheresignchange[0][ind]
			#pos = np.int_(xinterp[wheresignchange[0][ind]])
			#pos = xinterp[ind]
				
		#	pos = set[0][np.int_(xinterp[wheresignchange[0][0]])]
			print pos
		else :
			pos = 0
		
		position[i] = xinterp[pos]
		
		i+=1

	return position

def findFrontAmp_V1(data):
	amp = np.zeros(len(data))
	i=0
	for set in data:
		amp[i] = np.nanmax(set[1])
		i+=1
	return amp
	
def findFrontAmp(data):
	amp = np.zeros(len(data))
	i=0
	for set in data:
		m = [0,0]
		for j in range(2):
			m[j] = np.nanmax(set[1])
			set[1][np.nanargmax(set[1])] = 0
			
		amp[i] = np.nanmax(set[1])
		i+=1
	return amp

def frontSpeed(x,t):	
	return np.gradient(x[1:],np.diff(t))
	
def Cell2Meter(args):
	param = IO.Params(folder = args.folder).get()
	L = float(param["unit_l"])/3.085677e16
	dx = pow(2,- args.level)*L
	return dx
	
def plot_xft(args):
	data = readAllProfil(folder = args.folder)	
	t = gettmap(folder = args.folder)
	unit = Cell2Meter(args)
	y = findFrontPosition(data)
	y = np.multiply(y,unit)

	plt.plot(t[1:],y[1:])
#	plt.xlim(0,90)
#	plt.ylim(0,600)
	plt.ylabel(r'position (pc)')
	plt.xlabel(r'time (Myr)')
		
def plot_vft(args):
	data = readAllProfil(folder = args.folder)	
	t = gettmap(folder = args.folder)
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
	t = gettmap(folder = args.folder)
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
