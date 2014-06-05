#!/usr/bin/env python
import sys ,time
import numpy as np

import matplotlib.pylab as plt
from fonction_part import *
from fonction_IO import *
from fonction_physique import *
from evol import *

sys.path.append("/home/observatoire/deparis/Quartz/python/PyKDtree")
import kdtree as kd

def readSize(filename):
	f=open(filename,'r')
	npart=int(f.readline())
	npartingroup=int(f.readline())
	ngroup=int(f.readline())

	print 'npart={0} # npart in group={1} # ngroup={2}'.format(npart,npartingroup,ngroup)

	vsize=np.zeros(ngroup,dtype=np.int32)

	i=0
	for line in f:
		line=line.strip()
		column=line.split()
		vsize[i]=int(column[1])
		i=i+1

	f.closed
	return vsize	



def readDen(filename):

	denname = filename[:-10]  + "halo" + filename[-6:] + ".den"

	file = open(denname, "rb")
	print "Reading file", denname

	N 	= np.fromfile(file, dtype=np.int32   ,count=1)[0]
	den 	= np.fromfile(file, dtype=np.float32 ,count=N)
	
	file.close()

	print "Read OK"
	return N, den



def readTag(filename):

	tagname = filename[:-10]  + "halo" + filename[-6:] + ".tag"

	file = open(tagname, "rb")
	print "Reading file", tagname

	dummy 	= np.fromfile(file, dtype=np.int32   ,count=1)
	npart 	= np.fromfile(file, dtype=np.int32   ,count=1)[0]
	ngrp 	= np.fromfile(file, dtype=np.int32   ,count=1)[0]
	dummy 	= np.fromfile(file, dtype=np.int32   ,count=1)
	dummy 	= np.fromfile(file, dtype=np.int32   ,count=1)
	tag 	= np.fromfile(file, dtype=np.int32   ,count=npart)
	dummy 	= np.fromfile(file, dtype=np.int32   ,count=1)

	file.close()

	print "Read OK"
	return npart, ngrp, tag



def getCenter(halo, den):

	mask = np.argmax(den)
	x  = halo.x[mask]
	y  = halo.y[mask]
	z  = halo.z[mask]
	
	return np.array((x,y,z), dtype = np.float32)
	
def getRvir(halo):		

	M = np.sum(halo.mass)
	return np.power(3.*M/(800*np.pi),1./3.)

def getHalo(part, tag, den, HALO):

	mask  = np.where(tag == HALO)
	nmask = len(mask[0])

	halo = Part(nmask,0)
	halo.mask(part,mask)

	return halo, getCenter(halo, den[mask]), getRvir(halo)



def findStar(center,Rvir, star, tree):

	if tree : 
		mask3 = tree.locatenear( center, Rvir )	
		
	else : 
		mask3 = []

	nmask3 = len(mask3)
	starshalo = Part(nmask3,1)

	if nmask3 != 0:
		starshalo.mask(star,mask3)

	return starshalo, np.array(mask3, dtype = np.int32)

def genHopFiles(args):

	folder  = "../utils/hop/"


	for f in args.files:
		outname = f[:-10]  + "halo" + f[-6:]
	
		try :
			file = open(outname + ".den", "rb")
			file.close()
		except  IOError:
			commande="./" + folder + "hop -in " + f + " -o " + outname +" -p 1 -nf " + str(args.nproc)
			print commande
			os.system(commande)


		try :
			file = open(outname + ".tag", "rb")
			file.close()
		except  IOError:
			commande =  "./" + folder + "regroup -root output_hop -douter 80. -dsaddle 200. -dpeak 240. -f77 -o " + outname
			print commande
			os.system(commande)



def readAll(filename):
	
	
	npartden, den       = readDen(filename)
	nparttag, ngrp, tag = readTag(filename)
	npart,a,part        = readPart(filename  ,args)
	nstar,a,star        = readStars(filename ,args)
	
	return den, tag, part, star, nstar, a,ngrp

def genHaloFiles(args):

	for file in args.files:

		try :
			outname = file[:-10]  + "halo" + file[-6:] + ".halo"
			f = open(outname,'rb')	
			f.close()

		except  IOError:

			t0 = time.time()	
			den, tag, part, star, nstar, a, ngrp = readAll(file)	


			t1 = time.time()
			if nstar :	
				print "Tree generation"
				tree = kd.Tree(star.x,star.y,star.z)
				print "Tree generation OK"
			else :
				tree = 0


			t2 = time.time()
			print "Computation"


			center   = np.empty((ngrp), dtype=np.object)
			mask     = np.empty((ngrp), dtype=np.object)


			nstarhalo= np.zeros(ngrp, dtype=np.int32)
			Rvir     = np.zeros(ngrp, dtype=np.float32)
			Mh       = np.zeros(ngrp, dtype=np.float32)
			Ms       = np.zeros(ngrp, dtype=np.float32)


			for HALO in range(ngrp):

		
				halo,center[HALO],Rvir[HALO]   = getHalo(part,tag,den,HALO)
				starHalo, mask[HALO]           = findStar(center[HALO],Rvir[HALO],star, tree)

				nstarhalo[HALO] = len(mask[HALO])
				Mh[HALO] = np.sum(halo.mass)
				Ms[HALO] = np.sum(starHalo.mass)



			outname = file[:-10]  + "halo" + file[-6:] + ".halo"
			f = open(outname,'wb')		

			ngrp.tofile(f)		

			Rvir.tofile(f)
			Mh.tofile(f)
			Ms.tofile(f)
			nstarhalo.tofile(f)
			for i in range(ngrp):
				center[i].tofile(f)
				mask[i].tofile(f)			
			f.close()


			t3 = time.time()	

			print "Computation OK"
			print "Lecture", t1-t0
			print "tree",    t2-t1	
			print "Calcul",  t3-t2


def readHaloFile(file):


	outname = file[:-10]  + "halo" + file[-6:] + ".halo"
	try :
		f = open(outname,'rb')	
	except  IOError:
		print outname, "not found"

	ngrp = np.fromfile(f, dtype = np.int32, count = 1)		
	print ngrp

	center   = np.empty((ngrp), dtype=np.object)
	Rvir     = np.zeros( ngrp)
	Mh       = np.zeros( ngrp)
	Ms       = np.zeros( ngrp)	
	N       = np.zeros( ngrp)	
	mask     = np.empty((ngrp), dtype=np.object)


	Rvir   = np.fromfile(f, dtype = np.float32, count =ngrp)
	Mh     = np.fromfile(f, dtype = np.float32, count =ngrp)
	Ms     = np.fromfile(f, dtype = np.float32, count =ngrp)
	N      = np.fromfile(f, dtype = np.int32,   count =ngrp)

	for i in range(ngrp):
	
		center[i] = np.fromfile(f, dtype = np.float32, count = 3)
		mask[i]   = np.fromfile(f, dtype = np.int32, count =N[i])			

	f.close()


	return center,Rvir,Mh,Ms,N,mask




def SFR_f_Mh(Mh, mask, file) : 

	nstar,a,star        = readStars(file ,args)

	ngrp = len(mask)
	starshalo     = np.empty((ngrp), dtype=np.object)

	for grp in range(ngrp):

		nmask = len(mask[grp])
		starshalo[grp] = Part(nmask,1)
		if nmask != 0:


			starshalo[grp].mask(star,mask[grp])



	b = 1

	bin_edges = np.zeros(b+1)
	n= len(Mh)-1

	for i in range(0,b):
		bin_edges[i] = Mh[n/(i+1)]
	bin_edges[b] = Mh[0]

	print bin_edges


#	print	np.min(Mh)
#	n0, bin_edge = np.histogram(Mh, bins = bin_edges&)



	for i in range(b):

		mask0 = np.where( np.logical_and(   Mh > bin_edges[i], Mh < bin_edges[i+1]   ) )

		if len(mask0):	
			tmp = Part(0,1)			
			tmp.append(starshalo[mask0])

			SFRFromSnap(tmp, file + "\t"+ str(i))

	observation()

def Ms_f_Mh(Mh0,Ms0):

		L= 10
		Mh = m2mo(Mh0,L)
		Ms = m2mo(Ms0,L)


		b = 10   
		n0, bins0 = np.histogram(Mh,bins=b)
		n1, bins1 = np.histogram(Ms,bins=b)
		

		x = np.zeros(b)
		y = np.zeros(b)
		for i in range(b):
			x[i] = n0[i] * bins0[i] + ( bins0[i+1] - bins0[i] )/2.
			y[i] = n1[i] * bins1[i] + ( bins1[i+1] - bins1[i] )/2.

		print x,y

		plt.loglog(n0,n1, label=file)


		"""
		xE = np.zeros(b)
		for i in range(b):
			xE[i] =  bins0[i+1] - bins0[i]


		ax = plt.subplot(111)
		ax.set_xscale("log", nonposx='clip')
		ax.set_yscale("log", nonposy='clip')

		plt.errorbar(n0, n1 )
		"""
if __name__ == "__main__":

	args = getargs()

	genHopFiles(args)
	genHaloFiles(args)



	for file in args.files:
		center,Rvir,Mh,Ms,N,mask = readHaloFile(file)
	


		SFR_f_Mh(Mh, mask, file)
		#Ms_f_Mh(Mh,Ms)



		

#		x[i] = a2z(a)
#		y[i] = np.mean(Rap) + 1.
#		i+=1
		



#	plt.loglog(x,y,'.')
	plt.legend()
	plt.show()




	"""



		
	fig = plt.figure()	

	
		ax = fig.add_subplot(1, 1, 1)

		plt.plot(starHalo.x,starHalo.y, '.')
		
		ax.add_patch( plt.Circle((center[0],center[1]),Rvir,color='b',fill=False))
		plt.plot(center[0],center[1], 'ro')


	
	
	print out


	"""












