from fonction_part import *
from fonction_IO import *
import sys
import numpy as np
import matplotlib.pylab as plt

if __name__ == "__main__":

	args = getargs()
	foldername=args.folder
	files= args.files	
	
	
	N = np.zeros(len(files), dtype=np.int32)
	A = np.zeros(len(files), dtype=np.float64)
	Mtot = np.zeros(len(files), dtype=np.float64)

	for i in range(len(files)) :
		N[i] = getN(files[i])
		A[i] = getA(files[i]) 
		Mtot[i] = getMtotPart(files[i])

	print Mtot

	Fname = foldername + "t.dat"
	fileout = open(Fname, "w")
	print "writing file ", Fname

	fileout.write(	str(len(files)).zfill(8))

	N.tofile(fileout)
	A.tofile(fileout)
	Mtot.tofile(fileout)

	fileout.close()

