import sys
import os 
import numpy as np
import matplotlib.pylab as plt


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
	filePart = open(filename, "rb")
	return np.fromfile(filePart, dtype=np.int32, count=1)[0]

def getA(filename):
	filePart = open(filename, "rb")
	N = np.fromfile(filePart, dtype=np.int32, count=1)[0]
	t = np.fromfile(filePart, dtype=np.float32, count=1)[0]

	return t

def getZ(filename):
	return 1.0/getA(filename) -1

def getNproc(foldername):
	files = os.listdir(foldername)
	tmp =0
	for file in files:
		if file[0:10]=="part.00000" and file[-3:]!=".3D" :	
			tmp +=1
	return tmp

def ReadPart(filename):

	Ntot=0
	nProc = getNproc("../data")
#	nProc = 1
	parts = []

	for proc in range(nProc):
		filename = filename[0:-5] + str(proc).zfill(5)
		print "Reading file ", filename
		filePart = open(filename, "rb")
		N = np.fromfile(filePart, dtype=np.int32, count=1)[0]
		t = np.fromfile(filePart, dtype=np.float32, count=1)[0]
		Ntot += N
	
		for x in range(0,N) :
			parts.append(Part(filePart))

	return N,t,parts

def ReadStars(filename):

	N,t,parts=ReadPart(filename)
	stars=[]
	nstars=0

	for p in parts :
		if (p.isStar) :
			stars.append(p)
			nstars+=1


	print nstars
	return nstars,t,stars

#################################################################################################################################


def getStars(parts) : 
	stars=[]		
	n=0	
	for i in range (0,len(parts)):
		if parts[i].isStar :
			stars.append(parts[i])
			n+=1
	print n , "etoiles\n"
	return stars


def listPart(foldername):

	files = os.listdir(foldername)

	tmp=[]
	for file in files:
		if file[0:3]== "par" :	
			tmp.append(foldername + file)

	return  sorted(tmp)

def listOriginePart(foldername) :

	files = listPart(foldername)

	tmp=[]
	for file in files :
		if file[-3:]!= ".3D" and file[-6]== "p"  :
			tmp.append( file)
	return tmp

def makevideofile(foldername) :
	
	files = listPart(foldername)

	f = open("../data/stars.visit", "wt")

	tmp=[]
	for file in files:
		if file[0:3]== "par" :	
			tmp.append(file)

	for i in range(len(files)):
		if tmp[i][-9:]==".stars.3D" :
			f.write(files[i] )
			f.write( "\n" )
	f.close()	

	
def PartToVisit(parts, NameFileOut) :
	print "writing file", NameFileOut
	f = open(NameFileOut, "wt")
	f.write("x y z value\n");
	for p in parts:
	    f.write("%g %g %g %g\n" % (p.x,p.y,p.z,p.age))
	f.close()	

def FolderToVisit(foldername):

	files = listPart(foldername)

	for i in range(len(files)):
		file = files[i]
		if file[0:3]== "par" and file[-1]!="D":	
			N,t,parts= ReadPart(file)		
			PartToVisit(parts, file+".3D")
#			PartToVisit(getStars(parts), foldername + file +".stars.3D")

def getMtotPart(file):

	N,t,stars = ReadStars(file)
	M=0

	for s in stars :
		M += s.mass
	
	return M


