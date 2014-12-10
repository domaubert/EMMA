#!/usr/bin/env python
import sys, os
import numpy as np
import matplotlib.pylab as plt
from pylab import *
from numpy import ma
from operator import add


from fonction_IO  import *
from fonction_amr import *
from plotpart import *


if __name__ == "__main__":	


	folder = "data/movie/"
	files = os.listdir(folder)
	i,j=0,0
	
	nsub = 2
	fig = plt.figure()

	for file in np.sort(files):
		i+=1

		if  i%nsub == 0 :
			j+=1		
			print "%d / %d"%(j,len(files)/nsub)
			f = open(folder + file, "rb")
			x 	= np.fromfile(f, dtype=np.int32  ,count=1)[0]
			y 	= np.fromfile(f, dtype=np.int32  ,count=1)[0]
			a 	= np.fromfile(f, dtype=np.float32  ,count=1)[0]		

			m1 	= np.fromfile(f, dtype=np.float32  ,count=x*y)		#pot
			m2 	= np.fromfile(f, dtype=np.float32  ,count=x*y)		#den
			m3 	= np.fromfile(f, dtype=np.float32  ,count=x*y)		#X
			m4 	= np.fromfile(f, dtype=np.float32  ,count=x*y)		#temp

			data = m2.reshape(x,y)

			plt.imshow(data, interpolation='nearest')
			plt.colorbar()
		#	plt.clim(np.log(0.0135),np.log(0.0175))
			plt.savefig("data/img/"+file+".png")
			plt.clf()

			f.close()



