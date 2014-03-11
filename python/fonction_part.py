import sys
import os 
import numpy as np
import matplotlib.pylab as plt
from fonction_IO import *

class Part : 
	def __init__(self,filePart):
		self.x  = np.float64(np.fromfile(filePart, dtype=np.float32, count=1))
		self.y  = np.float64(np.fromfile(filePart, dtype=np.float32, count=1))
		self.z  = np.float64(np.fromfile(filePart, dtype=np.float32, count=1))
		self.vx = np.float64(np.fromfile(filePart, dtype=np.float32, count=1))
		self.vy = np.float64(np.fromfile(filePart, dtype=np.float32, count=1))
		self.vz = np.float64(np.fromfile(filePart, dtype=np.float32, count=1))	

		self.idx  = np.float64(np.fromfile(filePart, dtype=np.float32, count=1))
		self.mass = np.float64(np.fromfile(filePart, dtype=np.float32, count=1))
		self.epot = np.float64(np.fromfile(filePart, dtype=np.float32, count=1))
		self.ekin = np.float64(np.fromfile(filePart, dtype=np.float32, count=1))

		self.isStar = np.fromfile(filePart, dtype=np.int32, count=1)
		self.age    = np.float64(np.fromfile(filePart, dtype=np.float32, count=1))


	def printPos(self):	
		print self.x, self.y, self.z

	def printVit(self):	
		print self.vx, self.vy, self.vz


#################################################################################################################################

def getN(filename):
	try :
		filePart = open(filename, "rb")
	except IOError:
		print 'cannot open', filename 
		return 1
		
	return np.fromfile(filePart, dtype=np.int32, count=1)[0]
	filePart.close()


def getA(filename):
	try : 
		filePart = open(filename, "rb")
	except IOError:
			print 'cannot open', filename 

	N = np.fromfile(filePart, dtype=np.int32, count=1)[0]
	a = np.float64(np.fromfile(filePart, dtype=np.float32, count=1))
	filePart.close()
	return a

def getZ(filename):
	return 1.0/getA(filename) -1

def ReadPart(filename, args):

	Ntot=0
	nProc = getNproc(args)
	parts = []
	t=0

 
	print "Reading file ", filename[0:-7]
	for proc in range(nProc):
		filename = filename[0:-5] + str(proc).zfill(5)		

		try:
			filePart = open(filename, "rb")
		except IOError:
			print 'cannot open', filename 
	
		N = np.fromfile(filePart, dtype=np.int32, count=1)[0]
		t = np.fromfile(filePart, dtype=np.float32, count=1)[0]
		Ntot += N

		for x in range(0,N) :
			parts.append(Part(filePart))

		filePart.close()
	return Ntot,t,parts


	

def getStars(parts) : 
	stars=[]
	nstars=0

	for p in parts :
		if (p.isStar) :
			stars.append(p)
			nstars+=1
	print nstars, "stars"
	return nstars,stars


def ReadStars(filename, args):

	N,t,parts=ReadPart(filename,args)
	nstars, stars = getStars(parts)

	return nstars,t,stars


#################################################################################################################################

def getMtotPart(file, args):

	N,t,stars = ReadStars(file, args)
	M=0
	for s in stars :
		M += s.mass

	return M




