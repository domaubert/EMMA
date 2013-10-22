
import numpy as np
import matplotlib.pylab as plt

# ========================================================
# ========================================================
def read_size(filename):
    f=open(filename,'r')
    npart=int(f.readline())
    npartingroup=int(f.readline())
    ngroup=int(f.readline())

    print 'npart={0} # npart in group={1} # ngroup={2}'.format(npart,npartingroup,ngroup)

    vsize=np.zeros(ngroup,dtype=np.int32)

    i=0
    for line in f:
        line=line.strip()
        column=line.split()
        vsize[i]=int(column[1])
        i=i+1
    
    
    f.closed
    return vsize
    
# ========================================================
# ========================================================

def read_hmf(filename):
    f=open(filename,'r')
    hdr=f.readline()

    m=np.zeros(1000,dtype=np.float32)
    dndm=np.zeros(1000,dtype=np.float32)
    i=0
    for line in f:
        line=line.strip()
        column=line.split()
        m[i]=float(column[0]) 
        dndm[i]=float(column[2]) 
        i=i+1

    m=m[:i]
    dndm=dndm[:i]
    f.closed
    return m,dndm

# ========================================================
# ========================================================
if __name__ == "__main__":
    import sys

#    fib(int(sys.argv[1]))
    
    m,dndm=read_hmf("fmass/mf9.dat")
    vs=read_size("fmass/snap9.size")
    
    omegam=0.3
    h0=70.
    lbox=64.
    npart=256**3

    rhoc=3.*(h0*1e3/3.086e22)**2/8./np.pi/6.67e-11
    mpart=rhoc*(omegam)*(lbox/(h0/100.)*3.086e22)**3/2e30/npart*(h0/100)

    binm=np.logspace(8,15,128)
    h,b=plt.histogram(vs*mpart,binm)
    h=h/(lbox**3)/np.diff(binm)
    binok=(binm[:-1]+binm[1:])*0.5

    #plt.clf()
    plt.loglog(binok,h,marker='o',linestyle='None')
    plt.loglog(m,dndm)
    plt.show()
    
