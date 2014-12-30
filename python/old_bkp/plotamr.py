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



def plothisto(args,field):
	
	filename = args.files[0][:-7]
	filename = filename[:-10] +  "grid" +  filename[-6:] 
	denoct2grid(filename, args, field,0)

	a=array("cube")
#	data =  np.log10(a.getData())
	data =  a.getData()
	n, bins, patches = plt.hist(data, 100, normed=0, histtype='stepfilled', log=True)
	plt.show()


def plotdiffslice(args):

#	args.field = ["rfield.src"]
#	args.field = ["rfield.snfb"]
	args.field = ["rfield.temp"]
#	args.field = ["field.d"]


	denoct2grid(args.files[0], args,0)
	data1 = np.sum( cube(args.files[0] + ".cube").getData()  ,axis=0)


	denoct2grid(args.files[1], args,0)
	data2 = np.sum( cube(args.files[1] + ".cube").getData()  ,axis=0)

#	plt.grid(True)
	

 

	data = data1-data2
	data = np.log10(data+ 1) 
#	data = data + np.abs(np.min(data))
#	data = ((np.abs(data/np.mean(data) - 1.0 )))


	print np.sum(data)


	N = pow(2,args.level)
	plt.imshow( data, interpolation='nearest',extent=(0,N,0,N), origin='top' ,label =args.files[0] + " - " + args.files[1] )
	plt.colorbar()

	N,t,parts=readStars(args.files[0], args)
	plotpart(args,N,t,parts)

	plt.legend()
	plt.show()

def plotslice(args):

#	args.field = ["rfield.src"]
#	args.field = ["rfield.xion"]
#	args.field = ["rfield.temp"]
	args.field = ["field.d"]
	args.field = ["level"]

	filename = args.files[0][:-7]
	filename = args.files[0]
	filename = filename[:-10] +  "grid" +  filename[-6:] 

	print filename

	denoct2grid(filename, args,0)
#	data = np.log10(np.sum( cube(filename + ".cube").getData()  ,axis=0))
	data = cube(filename + ".cube").getData()


	N = pow(2,args.level)
	data = data[N/2,:,:]


	#plt.grid(True)

	#data=np.abs(data) 
	plt.imshow(data, interpolation='nearest',extent=(0,N,0,N), origin='top' )
	plt.colorbar()
#	plt.clim(5,7)
#	N,t,parts=readStars(args.files[0], args)
#	plotpart(args,N,parts)

	plt.show()


def plotv(args):

	filename = args.files[0][:-7]
	filename = filename[:-10] +  "grid" +  filename[-6:] 

	args.field = ["field.d"]
	denoct2grid(filename, args,0)
	data = cube(filename + ".cube").getData()

	data =  np.sum( data,axis=0)
	data = np.log10(data)

	plt.grid(True)
	
	N = pow(2,args.level)
	plt.imshow( data, interpolation='nearest',extent=(0,N,0,N), origin='top' )
	plt.colorbar()


 
	args.field = ["field.u"]
	denoct2grid(filename, args,0)
	vx = cube(filename + ".cube").getData()
	vx =  np.array(np.mean( vx,axis=0))
	

	args.field = ["field.v"]
	denoct2grid(filename, args,0)
	vy = cube(filename + ".cube").getData()
	vy =  np.array(np.mean( vy,axis=0))

	args.field = ["field.w"]
	denoct2grid(filename, args,0)
	vz = cube(filename + ".cube").getData()
	vz =  np.array(np.mean( vy,axis=0))


	X,Y = np.meshgrid( arange(0,128),arange(0,128 ))
	
	M = sqrt(pow(vx, 2) + pow(vy, 2) + pow(vz, 2) )*1000



        s = 2
        Q = plt.quiver(X[::s,::s], Y[::s,::s], vx[::s,::s], vy[::s,::s], M[::s,::s], pivot='mid')
	
	qk = plt.quiverkey(Q, 1000, 1.05, 1, r'$1 \frac{m}{s}$',
                            labelpos='E',
                            fontproperties={'weight': 'bold'})


	plt.show()


def plotTR(args):
	rho = getField(args,args.files[0], "field.d")
	T   = getField(args,args.files[0], "temp")

	plt.loglog(rho/np.mean(rho),T, '.')


def plotStarTR(args):

	data = readFile("RT")

	data = [x for x in data.split('\n')]

	for i in range(len(data)):
		data[i] = [x for x in data[i].split('\t')] 

	n    = np.zeros(len(data))
	rho  = np.zeros(len(data))
	temp = np.zeros(len(data))
	z = np.zeros(len(data))
	ns=0

	for i in range(len(data)-1):
		if ( len(data[i])==6 and isfloat(data[i][1]) and isfloat(data[i][2]) and isfloat(data[i][3])  ) :
			if float(data[i][1]):
				ns +=1
				n[i] = float(data[i][1])
				rho[i] = float(data[i][2])/0.154330706937
				temp[i] = float(data[i][3])
				z[i] = float(data[i][5])
	

	lim1 = 3
	lim2 = 2

	mask = np.where(z>lim1)
	plt.loglog(rho[mask], temp[mask], 'r.')
	mask = np.where(np.logical_and(z<lim1, z>lim2) )
	plt.loglog(rho[mask], temp[mask], 'y.')
	mask = np.where(z<lim2)
	plt.loglog(rho[mask], temp[mask], 'g.')



def plotdiag(args):

	plotTR(args)
#	plotStarTR(args)

#	plt.xlim(1e-2, 1e4)
#	plt.ylim(1e-2, 1e8)
	
	#plt.legend(getZ(args.files[0]))

	plt.show()




if __name__ == "__main__":	


	args = getargs()

	filename = args.files[0]
#	filename = filename[:-10] +  "part" +  filename[-6:] 

	print "redshift" , getZ(filename + ".p00000"  )


#	plothisto(args, field)
	plotslice(args)

#	plotdiffslice(args)
#	plotv(args)
#	plotdiag(args)








