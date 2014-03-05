import sys
import os 
import numpy as np
import matplotlib.pylab as plt


class Part : 
	def __init__(self,filePart):
		self.x  = np.fromfile(filePart, dtype=np.float32, count=1)
		self.y  = np.fromfile(filePart, dtype=np.float32, count=1)
		self.z  = np.fromfile(filePart, dtype=np.float32, count=1)
		self.vx = np.fromfile(filePart, dtype=np.float32, count=1)
		self.vy = np.fromfile(filePart, dtype=np.float32, count=1)
		self.vz = np.fromfile(filePart, dtype=np.float32, count=1)	

		self.idx  = np.fromfile(filePart, dtype=np.float32, count=1)	
		self.mass = np.fromfile(filePart, dtype=np.float32, count=1)
		self.epot = np.fromfile(filePart, dtype=np.float32, count=1)	
		self.ekin = np.fromfile(filePart, dtype=np.float32, count=1)

		self.isStar = np.fromfile(filePart, dtype=np.float32, count=1)		
		self.age    = np.fromfile(filePart, dtype=np.float32, count=1)

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



def ReadPart(filename):

	print "Reading file ", filename
	filePart = open(filename, "rb")
	N = np.fromfile(filePart, dtype=np.int32, count=1)[0]
	t = np.fromfile(filePart, dtype=np.float32, count=1)[0]
	
	parts = []
	for x in range(0,N) :
		parts.append(Part(filePart))

	return N,t,parts

def ReadStars(filename):

	print "Reading file ", filename
	filePart = open(filename, "rb")
	N = np.fromfile(filePart, dtype=np.int32, count=1)[0]
	t = np.fromfile(filePart, dtype=np.float32, count=1)[0]
	
	nstars = 0
	parts = []
	for x in range(0,N) :
		tmp = Part(filePart)
		if (tmp.isStar) :
			parts.append(tmp)
			nstars+=1
	print nstars , "stars"
	return N,t,parts

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
	    f.write("%g %g %g %d\n" % (p.x,p.y,p.z,p.age))
	f.close()	

def FolderToVisit(foldername):

	files = listPart(foldername)

	for i in range(len(files)):
		file = files[i]
		if file[0:3]== "par" and file[-1]!="D":	
			N,t,parts= ReadPart(file)		
			PartToVisit(parts, file+".3D")
#			PartToVisit(getStars(parts), foldername + file +".stars.3D")

def getMtot(file):

	N,t,stars = ReadStars(file)
	M=0

	for s in stars :
		M += s.mass[0]		
	
	return M


