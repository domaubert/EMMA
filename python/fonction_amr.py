import sys
import os 
import numpy as np
import matplotlib.pylab as plt


class array :
	def __init__(self,fileName):

		file = open(fileName, "rb")
		REAL = np.dtype('Float64')

		self.x  = np.fromfile(file, dtype=np.int32, count=1)
		self.y  = np.fromfile(file, dtype=np.int32, count=1)
		self.z  = np.fromfile(file, dtype=np.int32, count=1)
		self.a  = np.fromfile(file, dtype=np.float32, count=1)
		self.N  = self.x[0]*self.y[0]*self.z[0]
		self.data=np.zeros( self.N, dtype=REAL)

		for i in range(self.N) : 
			self.data[i]  = np.fromfile(file, dtype=REAL, count=1)
	#	print self.data.mean()
	#	print self.data.min()
				
	def getArray(self):
		return self.N, self.a,  self.data

	def getN(self):
		return self.N
	def geta(self):
		return self.a

class cube :
	def __init__(self,fileName):

		file = open(fileName, "rb")
		REAL = np.dtype('Float64')

		self.x  = np.fromfile(file, dtype=np.int32, count=1)
		self.y  = np.fromfile(file, dtype=np.int32, count=1)
		self.z  = np.fromfile(file, dtype=np.int32, count=1)
		self.a  = np.fromfile(file, dtype=np.float32, count=1)

		self.N  = (self.x[0]*self.y[0]*self.z[0])

		for k in range(0, self.z[0]) : 
			for j in range(0, self.y[0]) : 
				for i in range(0, self.x[0]) : 				
					self.data[i,j,k]  = np.fromfile(file, dtype=REAL, count=1)
					

	def getData(self):
		return self.data

def denoct2grid(n):

	o2g  = "./../utils/oct2grid "
	data = "../data/grid." + str(n).zfill(5)
	out  = "cube"

	os.system(o2g + data + " 6 101 " + out + " 1 0 -1 0 1 0 1 0 1 ")
#	os.system(o2g + data + " 6 1 " + out + " 1 0 -1 0 1 0 1 0 1 ")


def getMtotGaz(SnapNumber):


	denoct2grid(SnapNumber)


	den=array("cube")
	N, a, data = den.getArray()


	lmax=6
	vcell = pow(2,-3* lmax)
	
	for i in range(N):
		data[i]	 *= vcell 

	return data.sum()



