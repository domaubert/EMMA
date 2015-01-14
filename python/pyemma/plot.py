import matplotlib.pylab as plt

from convert import *
from amr import *

def getSlice(filename,level,force=0, nproc=IO.getNproc(), field="density", xmin=0, xmax=-1, ymin=0, ymax=-1, zmin=0, zmax=-1, log=False):
	N = pow(2,level)
	if xmax == -1 :
		xmax = N
	if ymax == -1 :
		ymax = N
	if zmax == -1 :
		zmax = N
		
	oct2grid(filename,level, force, nproc, field, xmin, xmax, ymin, ymax, zmin, zmax)

	data = cube(filename.replace("grid","cube"+ str(level)) + "." + field)

	print data.geta()
	data = data.getData()

	if log:
		data = np.log10(data)
	
	#data = np.max(data,axis=0)
	data = np.mean(data,axis=0)
	return data

def slice(filename,level,force=0, nproc=IO.getNproc(), field="density", xmin=0, xmax=-1, ymin=0, ymax=-1, zmin=0, zmax=-1, log=False):
	data = getSlice(filename,level,force, nproc, field, xmin, xmax, ymin, ymax, zmin, zmax, log)

	plt.clf()
	plt.imshow(data, interpolation='nearest',extent=(xmin,xmax,ymin,ymax), origin='top' )
	plt.colorbar()
	plt.show(block=False)

def diffslice(filename1,filename2,level,force=0, nproc=IO.getNproc(), field="density", xmin=0, xmax=-1, ymin=0, ymax=-1, zmin=0, zmax=-1, log=False):
	data1 = getSlice(filename1,level,force, nproc, field, xmin, xmax, ymin, ymax, zmin, zmax, log)
	data2 = getSlice(filename2,level,force, nproc, field, xmin, xmax, ymin, ymax, zmin, zmax, log)

	data=np.subtract(data1,data2)

	"""
	if log:
		data = np.log10(np.abs(data))
	"""
	plt.clf()
	plt.imshow(data, interpolation='nearest',extent=(xmin,xmax,ymin,ymax), origin='top' )
	plt.colorbar()
	plt.show(block=False)
