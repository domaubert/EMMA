import sys, os
from struct import *
import numpy as np
import matplotlib.pylab as plt

class array :
	def __init__(self,fileName):
		fileName = fileName[:-15] +  "grid" +  fileName[-11:] 
		file = open(fileName, "rb")

		self.x,self.y,self.z    = np.fromfile(file, dtype=np.int32    ,count=3)
		self.a    		= np.fromfile(file, dtype=np.float32  ,count=1)[0]

		self.Z    = 1.0/self.a-1.0
		self.N    = int(self.x*self.y*self.z)
	
		self.data    		= np.fromfile(file, dtype=np.float32  ,count=self.N)

		
	def getData(self):
		return self.data
	def getN(self):
		return self.N
	def getSize(self):
		return self.x,self.y,self.z
	def geta(self):
		return self.a
	def getz(self):
		return self.z

class cube :
	def __init__(self,fileName):

		self.a    = array(fileName)	
		self.data = np.reshape(self.a.getData(),  (self.a.getSize()) ) 
		
	def getData(self):
		return self.data
	def getZ(self):
		return self.a.getZ()

def denoct2grid(data,args,silo):

	o2g  = "utils/oct2grid "

	out  = data + ".cube "
	f= "0 "
	if args.field[0]=="density" : 
		f=" 1 "
	if args.field[0]=="pot" : 
		f=" 2 "
	if args.field[0]=="cpu" : 
		f=" 3 "
	if args.field[0]=="marked" : 
		f=" 4 "
	if args.field[0]=="field.d" : 
		f=" 101 "
	if args.field[0]=="field.u" : 
		f=" 102 "
	if args.field[0]=="field.v" : 
		f=" 103 "
	if args.field[0]=="field.w" : 
		f=" 104 "
	if args.field[0]=="field.p" : 
		f=" 105 "
	if args.field[0]=="field.E" : 
		f=" 106 "
	if args.field[0]=="rfield.src" : 
		f=" 705 "
	if args.field[0]=="rfield.temp" : 
		f=" 707 "
	if args.field[0]=="field.u" : 
		f=" 102 "
#	if args.field[0]=="rfield.snfb" : 
#		f=" 706 "
	if args.field[0]=="rfield.xion" : 
		f=" 706 "
	if args.field[0]=="f0" : 
		f=" 201 "
	

	if f=="0 ":
		print "entrez un champ"
		sys.exit()

	print "Lecture de la grille", data 
	commande = o2g + " " + data + " " +  str(args.level) + " "  + str(f) + out + str(args.nproc) + " " + str(silo) + " -1 0 1 0 1 0 1 > dump" 
	print commande
	os.system(commande)


def getField(args, filename, field):
	filename = filename[:-7]
	filename = filename[:-10] +  "grid" +  filename[-6:] 

	args.field = [field]

	denoct2grid(filename, args,0)
	return array(filename + ".cube").getData()

def genCube(args, field):
	for file in args.files:
		filename = file[:-7]
		filename = filename[:-10] +  "grid" +  filename[-6:] 

		args.field = [field]
		denoct2grid(filename, args,0)


def getMtotGaz(folder,SnapNumber):

	denoct2grid(folder, SnapNumber, "field.d")
	den=array("cube")
	N, a, data = den.getArray()

	lmax=6
	vcell = pow(2,-3* lmax)
	
	for i in range(N):
		data[i]	 *= vcell 

	return data.sum()



