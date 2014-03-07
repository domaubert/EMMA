import sys
import os
import numpy as np
import matplotlib.pylab as plt

from fonction_part import *
from fonction_amr import *



if __name__ == "__main__":	


	if len(sys.argv)==2 :
		SnapNumber = sys.argv[1]

		Mtotgaz = getMtotGaz(SnapNumber)
		MtotPart =  getMtotPart("../data/part."  + str(SnapNumber).zfill(5) + ".p00000")

		print "================================ "
		print "Total mass of baryon wanted\t0.154330706937" 
		print "Total mass of baryon \t\t",  MtotPart + Mtotgaz
		print "Total mass of gaz \t\t",  Mtotgaz
		print "Total mass of Stars \t\t",  MtotPart
		print "================================ \n"



	elif len(sys.argv)==3 :
	
		n = np.arange(int(sys.argv[1]), int(sys.argv[2]))
		Mdata = np.zeros(len(n))

		for i in n :			
			Mdata[i-n[0]] = getMtotGaz(i) + getMtotPart("../data/part."  + str(i).zfill(5) + ".p00000")

		plt.plot(n,Mdata)
		plt.show()
