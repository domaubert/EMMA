#!/usr/bin/env python
import time
import sys
import numpy as np

from fonction_part import *
from fonction_IO import *
from decimal import Decimal


if __name__ == "__main__":

	args = getargs()
	foldername=args.folder
	files= args.files	

	n = 1

	size=  len(files)/n

	print size
	N = np.zeros(size, dtype=np.int32)
	A = np.zeros(size, dtype=np.float64)
	Mtot = np.zeros(size, dtype=np.float64)
	

	for j in range(0, len(files), n) :
		file = files[j][:-17] +  "star" +  files[j][-13:]
		N[j/n] = getNtot(file,args)
		A[j/n] = getA(file) 
		Mtot[j/n] = getMtotPart(file, args)


	print Mtot

	Fname = foldername[0] + "tdata.00000.p00000"
	fileout = open(Fname, "w")
	print "writing file ", Fname

	
	fileout.write(	str(size).zfill(8))

	N.tofile(fileout)
	A.tofile(fileout)
	Mtot.tofile(fileout)

	fileout.close()

