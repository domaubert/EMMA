from fonction_part import *
from fonction_IO import *
from fonction_physique import *
import numpy as np




if __name__ == "__main__":


	args=getargs()

	fs = open(args.folder + "t.dat")
	
	n = (int)(fs.read(8))
	N = np.fromfile(fs, dtype=np.int32, count=n) 
	A = np.fromfile(fs, dtype=np.float64, count=n) 
	Mtot = np.fromfile(fs, dtype=np.float64, count=n) 
	fs.close()

	npart = N[0]
	for i in range(n) :
		N[i] -= npart

	Z = np.zeros(n, dtype = np.float64)
	T =  np.zeros(n, dtype = np.float64)
	SFR = np.zeros(len(Mtot))

	L = 10 #Mpc

	for i in range(n) :
		Z[i] = a2z(A[i])
		T[i] = a2t(A[i]) 





	for i in range(n-1) :
		SFR[i] = m2mo(Mtot[i+1] - Mtot[i], A[i],L) / (T[i+1] - T[i]) / pow(L,3)




	f= findfirststar(N)+1
	print f

	print Mtot[f:].min(), m2mo(Mtot[f:],A[f:], L).min()
	print Mtot[f:].max(), m2mo(Mtot[f:],A[f:], L).max()



	plt.semilogy(Z[f:-1],SFR[f:-1])
	plt.xlabel(r'$Z$')
	plt.ylabel(r'$SFR (M_{o}.yr^{-1}.Mpc^{-3})$' )


	plt.show()











