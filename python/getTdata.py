#!/usr/bin/env python
import time
import sys
import numpy as np

from fonction_part import *
from fonction_IO import *


if __name__ == "__main__":

	args = getargs()
	foldername=args.folder
	files= args.files	


	
	n = 1
#	files=files[0:1]
	size = len(files)/n
	
	N = np.zeros(size, dtype=np.int32)
	A = np.zeros(size, dtype=np.float64)
	Mtot = np.zeros(size, dtype=np.float64)
	

	i=0
	for j in range(size) :
		file = files[i][:-17] +  "star" +  files[i][-13:]
		N[j] = getNtot(file,args)
		A[j] = getA(file) 
		Mtot[j] = getMtotPart(file, args)
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

