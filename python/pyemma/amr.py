import numpy as np

class array :
	def __init__(self,fileName):
		file = open(fileName, "rb")

		self.x,self.y,self.z = np.fromfile(file, dtype=np.int32    ,count=3)
		self.a = np.fromfile(file, dtype=np.float32  ,count=1)[0]

		self.Z = 1.0/self.a-1.0
		self.N = int(self.x*self.y*self.z)
	
		self.data = np.fromfile(file, dtype=np.float32  ,count=self.N)
		file.close()
		
	def getData(self):
		return self.data
	def getN(self):
		return self.N
	def getSize(self):
		return self.z,self.y,self.x
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
	def geta(self):
		return self.a.geta()
	def getSize(self):
		return self.a.getSize()
"""
def geta(fileName):
	file = open(fileName, "rb")
	x,y,z = np.fromfile(file, dtype=np.int32    ,count=3)
	a = np.fromfile(file, dtype=np.float32  ,count=1)[0]
	file.close()
	return a
"""
def geta(fileName):
	file = open(fileName, "rb")
	a = np.fromfile(file, dtype=np.float64  ,count=1)[0]
	file.close()
	return a
	
def cubename(filename,level,field):
	return filename.replace("grid","cube"+ str(level)) + "." + str(field)


