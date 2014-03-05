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
				
	def getData(self):
		return self.data
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


