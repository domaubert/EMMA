import sys
import os
import numpy as np
import matplotlib.pylab as plt

from fonction_part import *
from fonction_amr import *






def getMtotBaryon(SnapNumber) :

	Mtotgaz = getMtotGaz(SnapNumber)
	print "Gaz " , Mtotgaz
	MtotPart =  getMtotPart("../data/part."  + str(SnapNumber).zfill(5) + ".p00000")
	print "Parts " , MtotPart

	return MtotPart + Mtotgaz


if __name__ == "__main__":	


	if len(sys.argv)==2 :
		SnapNumber = sys.argv[1]
		M = getMtotBaryon(SnapNumber)
#		M = getMtotGaz(SnapNumber)
		print "\n ================================ \n Total mass",  M , "\n ================================ \n"

	elif len(sys.argv)==3 :
	
		n = np.arange(int(sys.argv[1]), int(sys.argv[2]))
		Mdata = np.zeros(len(n))
		Mdata2 = np.zeros(len(n))


		for i in range(len(n)) :
			Mdata[i] = getMtotBaryon(n[i])
#		for i in range(len(Mdata) - 1) :
#			Mdata2[i] = Mdata[i+1] - Mdata[i]

		plt.plot(n,Mdata)
		plt.show()
