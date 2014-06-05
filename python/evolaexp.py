#!/usr/bin/env python
from fonction_part import *
from fonction_IO import *
from fonction_physique import *
import numpy as np
import matplotlib.pylab as plt 
from matplotlib.ticker import NullFormatter


def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False

if __name__ == "__main__":


	filename = "out"
	try:
		f = open(filename, "r")	
	except IOError:
		print 'cannot open', filename 

	data = f.read()
	f.close()


	data = [x for x in data.split('\n')]


	for i in range(len(data)):
		data[i] = [x for x in data[i].split(' ')] 
	print len(data)
		
	y= np.zeros(len(data)-2)
	x= range(len(y))
	
	for i in range(len(x)): 
		y[i] =  float(data[i+1][3][5:]) - float(data[i][3][5:])
	


	plt.semilogy(x, y)


	plt.show()



