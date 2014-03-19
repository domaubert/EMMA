#!/usr/bin/env python
from fonction_part import *
from fonction_IO import *
import sys
import numpy as np
import matplotlib.pylab as plt

if __name__ == "__main__":

	args = getargs()
	foldername=args.folder
	files= args.files	

	n = 1
	size = len(files)/n 
	
	N = np.zeros(size, dtype=np.int32)
	A = np.zeros(size, dtype=np.float64)
	Mtot = np.zeros(size, dtype=np.float64)
	

	i=0
	for j in range(size) :
		N[j] = getNtot(files[i],args)
		A[j] = getA(files[i]) 
		Mtot[j] = getMtotPart(files[i], args)
		i+=n

	print Mtot

	Fname = foldername[0] + "tdata.00000.p00000"
	fileout = open(Fname, "w")
	print "writing file ", Fname

	fileout.write(	str(size).zfill(8))

	N.tofile(fileout)
	A.tofile(fileout)
	Mtot.tofile(fileout)

	fileout.close()

