import sys, os
from struct import *
import numpy as np
import matplotlib.pylab as plt


class array :
	def __init__(self,fileName):

		file = open(fileName, "rb")

		self.x    = np.fromfile(file, dtype=np.int32, count=1)
		self.y    = np.fromfile(file, dtype=np.int32, count=1)
		self.z    = np.fromfile(file, dtype=np.int32, count=1)
		self.a    = np.fromfile(file, dtype=np.float32, count=1)
		self.N    = self.x[0]*self.y[0]*self.z[0]

		self.data = np.zeros( self.N, dtype=np.float64)
		for i in range(self.N) : 
			self.data[i]  = np.fromfile(file, dtype=np.float64, count=1)

	#	print self.N
	#	print self.data.mean()
	#	print self.data.min()
				
	def getData(self):
		return self.data
	def getN(self):
		return self.N
	def geta(self):
		return self.a

class cube :
	def __init__(self,fileName):

		file = open(fileName, "rb")


		self.x  = np.fromfile(file, dtype=np.int32, count=1)
		self.y  = np.fromfile(file, dtype=np.int32, count=1)
		self.z  = np.fromfile(file, dtype=np.int32, count=1)
		self.a  = np.fromfile(file, dtype=np.float32, count=1)

		self.N  = self.x[0]*self.y[0]*self.z[0]

		self.data = np.zeros( (self.x,self.y,self.z), dtype=np.float32)

		for k in range(0, self.z[0]) : 
			for j in range(0, self.y[0]) : 
				for i in range(0, self.x[0]) : 				
					self.data[i,j,k]  =  np.fromfile(file, dtype=np.float64, count=1)
					
#np.fromfile(file, dtype=np.float32, count=1)			

	def getData(self):
		return self.data


def denoct2grid(data,args,silo):

	o2g  = "utils/oct2grid "

	out  = " cube "
	f= "0 "
	if args.field[0]=="density" : 
		f=" 1 "
	if args.field[0]=="pot" : 
		f=" 2 "
	if args.field[0]=="cpu" : 
		f=" 3 "
	if args.field[0]=="marked" : 
		f=" 4 "
	if args.field[0]=="temp" : 
		f=" 707 "
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
	if f=="0 ":
		print "entrez un champ"
		sys.exit()

	print "Lecture de la grille", data 
	commande = o2g + " " + data + " " +  str(args.level) + " "  + str(f) + out + str(args.nproc) + " " + str(silo) + " -1 0 1 0 1 0 1 > dump" 
	print commande
	os.system(commande)



def getMtotGaz(folder,SnapNumber):


	denoct2grid(folder, SnapNumber, "field.d")
	den=array("cube")
	N, a, data = den.getArray()

	lmax=6
	vcell = pow(2,-3* lmax)
	
	for i in range(N):
		data[i]	 *= vcell 

	return data.sum()



