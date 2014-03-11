import sys
import os
import numpy as np
import matplotlib.pylab as plt

from fonction_part import *
from fonction_amr import *
from fonction_IO import *


if __name__ == "__main__":	

	args = getargs()
	folder = args.folder

	if not (args.plot) :
		for i in range(len(args.files)):



			Mtotgaz = getMtotGaz(folder,args.snap[i])
			MtotPart =  getMtotPart(args.files[i], args)

			print "================================ "
			print "Total baryonic mass wanted\t0.154330706937" 
			print "Total barionic mass  \t\t",  MtotPart + Mtotgaz
			print "Delta \t\t\t\t", MtotPart + Mtotgaz -0.154330706937
			print "================================="
			print "Total mass of gaz \t\t",  Mtotgaz
			print "Total mass of Stars \t\t",  MtotPart
			print "================================ \n"

	else  :
		
		n=args.snap
		Mdata = np.zeros(len(n))

		print n, Mdata

		for i in n :			
			Mdata[i-n[0]] = getMtotGaz(folder, i) + getMtotPart(folder + "part."  + str(i).zfill(5) + ".p00000", args)

		plt.plot(n,Mdata)
		plt.show()
