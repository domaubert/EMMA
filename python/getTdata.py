from fonction_part import *
import numpy as np
import matplotlib.pylab as plt

if __name__ == "__main__":

	
	if len(sys.argv)>1 :
		foldername = sys.argv[1]
	else :
		foldername="../data/"
	

	

	files = listOriginePart(foldername)
	files = files[10:40]

	
	N = np.zeros(len(files), dtype=np.int32)
	A = np.zeros(len(files), dtype=np.float32)
	Mtot = np.zeros(len(files), dtype=np.float32)

	for i in range(len(files)) :
		N[i] = getN(files[i])
		A[i] = getA(files[i]) 
		Mtot[i] = getMtotPart(files[i])



	Fname = foldername + "t.dat"
	fileout = open(Fname, "w")
	print "writing file ", Fname

	fileout.write(	str(len(files)).zfill(8))

	N.tofile(fileout)
	A.tofile(fileout)
	Mtot.tofile(fileout)

	fileout.close()

