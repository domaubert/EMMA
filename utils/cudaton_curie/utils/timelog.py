#!/usr/bin/env python

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io.numpyio import fwrite, fread

def main():

    
    fichier=sys.argv[1]
    ttab=np.loadtxt(fichier)

    print ttab.shape
    print ttab[9000,:]
    plt.yscale('log')
    plt.plot(ttab[:,0],ttab[:,1],'k--',label='transport')
#    plt.plot(ttab[:,0],ttab[:,2],'r',label='chem')
    plt.plot(ttab[:,0],ttab[:,3],'b',label='ion+cool')
 #   plt.plot(ttab[:,0],ttab[:,4],'g',label='update')
    plt.plot(ttab[:,0],ttab[:,5],'r',label='bnd')
    plt.plot(ttab[:,0],ttab[:,6],'g',label='io')
    plt.plot(ttab[:,0],ttab[:,7],'k',label='total')
    plt.legend(loc=4)
    plt.show()
if __name__ == '__main__':
	main()
        
