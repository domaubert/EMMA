#!/usr/bin/env python

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io.numpyio import fwrite, fread

def main():

    NBND=16;
    NBND2=2*NBND
    
    dirout=sys.argv[1]
    NCELL=int(sys.argv[2])

    NGPUX=int(sys.argv[3])
    NGPUY=int(sys.argv[4])
    NGPUZ=int(sys.argv[5])

    ISNAP=int(sys.argv[6])

    NGPUTOT=NGPUX*NGPUY*NGPUZ

    NX=NCELL/NGPUX
    NY=NCELL/NGPUY
    NZ=NCELL/NGPUZ


    x=np.zeros([NCELL,NCELL,NCELL])
    t=np.zeros([NCELL,NCELL,NCELL])

    for k in range(NGPUZ):
        for j in range(NGPUY):
            for i in range(NGPUX):
                idx=i+j*NGPUX+k*NGPUX*NGPUY

                fname=dirout+'/snap.%05d.p%05d'%(ISNAP,idx)
                print 'reading '+fname

                fbuf=open(fname,mode="rb")
                nc=fread(fbuf,1,'i')
                tt=fread(fbuf,1,'f')
                xtemp=fread(fbuf,(NX+NBND2)*(NY+NBND2)*(NZ+NBND2),'f')
                ttemp=fread(fbuf,(NX+NBND2)*(NY+NBND2)*(NZ+NBND2),'f')
                fbuf.close()

                xtemp=xtemp.reshape((NZ+NBND2,NY+NBND2,NX+NBND2))
                ttemp=ttemp.reshape((NZ+NBND2,NY+NBND2,NX+NBND2))

                
                x[k*NZ:(k+1)*NZ,j*NY:(j+1)*NY,i*NX:(i+1)*NX]=xtemp[NBND:NBND+NZ,NBND:NBND+NY,NBND:NBND+NX]
                t[k*NZ:(k+1)*NZ,j*NY:(j+1)*NY,i*NX:(i+1)*NX]=ttemp[NBND:NBND+NZ,NBND:NBND+NY,NBND:NBND+NX]


    xx=np.linspace(-6.6,6.6,NCELL)
    plt.title('time='+str(tt))
    plt.subplot(221)
    plt.imshow(1.-x[NCELL/2,:,:])
    plt.subplot(222)
    plt.plot(xx,1.-x[NCELL/2,NCELL/2,:])
    plt.xlim(-6.6,6.6)
    plt.yscale('log')

    plt.subplot(223)
    plt.imshow(t[NCELL/2,:,:])
    plt.subplot(224)
    plt.plot(xx,t[NCELL/2,NCELL/2,:])
    plt.yscale('log')
    plt.xlim(-6.6,6.6)
    plt.ylim(5000.,3e4)
    plt.show()
                
if __name__ == '__main__':
	main()
        
