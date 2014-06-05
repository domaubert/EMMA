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


	filename = "RT"
	try:
		f = open(filename, "r")	
	except IOError:
		print 'cannot open', filename 

	data = f.read()
	f.close()


	data = [x for x in data.split('\n')]
	print len(data)

	for i in range(len(data)):
		data[i] = [x for x in data[i].split('\t')] 

	n    = np.zeros(len(data))
	rho  = np.zeros(len(data))
	temp = np.zeros(len(data))
	ns=0

	for i in range(len(data)-1):
		if ( len(data[i])==4 and isfloat(data[i][1]) and isfloat(data[i][2]) and isfloat(data[i][3])  ) :
			if float(data[i][1]):
				ns +=1
				n[i] = float(data[i][1])
				rho[i] = float(data[i][2])/ 0.154330706937
				temp[i] = float(data[i][3])

	print ns
	plt.loglog(rho, temp, '.')

	'''
	nullfmt   = NullFormatter()         # no labels

	# definitions for the axes
	left, width = 0.1, 0.65
	bottom, height = 0.1, 0.65
	bottom_h = left_h = left+width+0.02

	rect_scatter = [left, bottom, width, height]
	rect_histx = [left, bottom_h, width, 0.2]
	rect_histy = [left_h, bottom, 0.2, height]

	# start with a rectangular Figure
	plt.figure(1, figsize=(8,8))

	axScatter = plt.axes(rect_scatter)
	axHistx = plt.axes(rect_histx)
	axHisty = plt.axes(rect_histy)

	# no labels
	axHistx.xaxis.set_major_formatter(nullfmt)
	axHisty.yaxis.set_major_formatter(nullfmt)

	# the scatter plot:
	axScatter.scatter(x, y)

	# now determine nice limits by hand:
	binwidth = 0.25
	xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )
	lim = ( int(xymax/binwidth) + 1) * binwidth

	axScatter.set_xlim( (-lim, lim) )
	axScatter.set_ylim( (-lim, lim) )

	bins = np.arange(-lim, lim + binwidth, binwidth)
	axHistx.hist(x, bins=bins)
	axHisty.hist(y, bins=bins, orientation='horizontal')

	axHistx.set_xlim( axScatter.get_xlim() )
	axHisty.set_ylim( axScatter.get_ylim() )
	'''
	plt.show()



