#!/usr/bin/env python
import sys, os
import numpy as np
import matplotlib.pylab as plt
from fonction_IO import *

def readE():
	data = readFile("emissivity.out")

	data = [x for x in data.split('\n')]

	for i in range(len(data)):
		data[i] = [x for x in data[i].split(' ')] 
	
	z = data[20]
	z = [float(x) for x in z[1:]] 

	data = data[21:-1]
	for i in range(len(data)):		
		data[i] = [float(x) for x in data[i][1:]] 


	return np.array(z), np.array(data)



if __name__ == "__main__":

	z,data = readE()

	n= len(data)
	x=np.zeros(n)
	y=np.zeros(n)		

	m=np.zeros(60)		

	for t in range(60):	
		for i in range(n):
			x[i] = data[i][0]
			y[i] = data[i][t]


	#	m[t] = y[i].mean()

		plt.semilogy(z, y[t])

	plt.show()
