
import numpy as np
import matplotlib.pylab as plt
import hop as hop
from mpl_toolkits.mplot3d import Axes3D

# ========================================================
# ========================================================
# ========================================================

class Part:

#==================================================    
    def loadfile_mono(self,fname):
        print 'PHASE LOAD reading '+fname
        f=open(fname,"rb")
        self.npart=np.fromfile(f,dtype=np.int32,count=1)[0]
        self.time=np.fromfile(f,dtype=np.float32,count=1)[0]
        self.x=np.zeros(self.npart,dtype=np.float32)
        self.y=np.zeros(self.npart,dtype=np.float32)
        self.z=np.zeros(self.npart,dtype=np.float32)
        self.vx=np.zeros(self.npart,dtype=np.float32)
        self.vy=np.zeros(self.npart,dtype=np.float32)
        self.vz=np.zeros(self.npart,dtype=np.float32)
        self.id=np.zeros(self.npart,dtype=np.float32)
        

        temp=np.fromfile(f,dtype=np.float32,count=7*self.npart)
        temp=(temp.reshape([self.npart,7])).transpose()
        self.x=temp[0,]
        self.y=temp[1,]
        self.z=temp[2,]
        self.vx=temp[3,]
        self.vy=temp[4,]
        self.vz=temp[5,]
        self.id=temp[6,]
        f.close

#==================================================    
    def getnpart(self,fname):
        print 'NPART LOAD reading '+fname
        f=open(fname,"rb")
        self.npart=np.fromfile(f,dtype=np.int32,count=1)[0]
        f.close


#==================================================    
    def alloc(self,npart):
        self.x=np.zeros(npart,dtype=np.float32)
        self.y=np.zeros(npart,dtype=np.float32)
        self.z=np.zeros(npart,dtype=np.float32)
        self.vx=np.zeros(npart,dtype=np.float32)
        self.vy=np.zeros(npart,dtype=np.float32)
        self.vz=np.zeros(npart,dtype=np.float32)
        self.id=np.zeros(npart,dtype=np.float32)
        self.npart=npart
        self.time=0.
#==================================================    
    def loadfile(self,fname,ncpu):
        
        nparttot=0
        icpu=0
        while icpu<ncpu:
            p=Part()
            p.getnpart(fname+'.p{0:05d}'.format(icpu))
            nparttot=nparttot+p.npart
            icpu=icpu+1

        print '{0} particles found in {1} file'.format(nparttot,ncpu)
        self.alloc(nparttot)
        icpu=0
        ncur=0

        while icpu<ncpu:
            p=Part()
            p.loadfile_mono(fname+'.p{0:05d}'.format(icpu))


            self.x[ncur:ncur+p.npart]=p.x
            self.y[ncur:ncur+p.npart]=p.y
            self.z[ncur:ncur+p.npart]=p.z

            self.vx[ncur:ncur+p.npart]=p.vx
            self.vy[ncur:ncur+p.npart]=p.vy
            self.vz[ncur:ncur+p.npart]=p.vz

            ncur=ncur+p.npart
            icpu=icpu+1

        
#==================================================    

    def loadtag(self,fname):
        print 'TAG LOAD reading '+fname
        f=open(fname,"rb")
        dummy=np.fromfile(f,dtype=np.int32,count=5)
        self.tag=np.fromfile(f,dtype=np.int32,count=self.npart)
        f.close

#==================================================    
    def __init__(self):
        self.x = []
        self.y = []
        self.z = []
        self.vx = []
        self.vy = []
        self.vz = []
        self.id = []
        self.npart=[]
        self.time=[]
        self.tag=[]

# ========================================================
# ========================================================
# ========================================================


if __name__ == "__main__":
    import sys
    
    
    i=0
    p=Part()
    p.loadfile('/home/daubert/Project/Quartz/data/part.00023',4)
    p.loadtag('/home/daubert/Project/Quartz/utils/hop/fmass/snap23.tag')

    w=np.where(p.tag==100)
    d,c,b=plt.histogram2d(p.x[w],p.y[w],bins=256,range=[[p.x[w].min(),p.x[w].max()],[p.y[w].min(),p.y[w].max()]])

    
    plt.clf()
    plt.imshow(np.log(1+d))
    plt.show()
    

    '''
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(p.x[w], p.y[w], p.z[w])
    '''

    plt.show()
