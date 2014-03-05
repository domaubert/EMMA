import sys
import os
import numpy as np
import matplotlib.pylab as plt
from fonction_part import *



def PlotHisto2d(parts):
	N=len(parts)
	x=np.zeros(N)
	y=np.zeros(N)
	for i in range(N) :
		x[i]=parts[i].x
		y[i]=parts[i].y

	Lmax=6
	H, xedges, yedges = np.histogram2d(x, y, bins=(pow(2,Lmax),pow(2,Lmax)))


	plt.imshow(H, interpolation='nearest')
	plt.colorbar()
	plt.show()


def plothisto2Dstars() :
	if len(sys.argv)>1 :
		filename = "../data/part." + sys.argv[1].zfill(5) + ".p00000"
		N,t,parts= ReadPart(filename)

		PlotHisto2d( getStars(parts) )


	else :
		print "Quel plot???"



def spectre() :

	if len(sys.argv)>1 :
		filename = "../data/part." + sys.argv[1].zfill(5) + ".p00000"
		N,t,parts= ReadPart(filename)

		stars=[]

		for i in range (0,N):
			if parts[i].isStar :
				stars.append(parts[i].mass[0])


		plt.hist(stars,50 , normed=1, log=True,  histtype='stepfilled')#,  cumulative=True


	else :
		print "Quel plot???"	

def plotpart() :
	if len(sys.argv)>1 :
		filename = "../data/part." + sys.argv[1].zfill(5) + ".p00000"

		N,t, stars=ReadStars(filename)
		

		plt.hist(stars, 50 , normed=1, histtype='stepfilled')#, log=False,  cumulative=True


	else :
		print "Quel plot???"	
	


if __name__ == "__main__":
	



#	plotpart()
#	plothisto2Dstars()
	spectre()	

	



	plt.show()
