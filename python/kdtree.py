"""
Module to create and use KD-trees: 
all you ever ever need for your neighbours searches.

How to make it work:

This module comes with at least
- This file ( duh. )
- The shared c library libkdtree.so

The quick way is to store the libkdtree.so somewhere, 
update the path at the line with ctypeslib.load_library, 
and voila. ( hopefully )



The hard way if you wish to modify the core c code:

The complete package should contain:
- The python file kdtree.py
- A readme file
- the makefile
- two libraries: libkdtree.so and libsort.so
- two directories: argsort and kdtree
     - argsort:
         argsort.c
         argsort.h
     - kdtree:
         kdtree.c
         kdtree.h


  The core functions building and using the kd tree are in
kdtree.c. These use an argsort function defined in
argsort.c . So if you want to build the libkdtree.so, you
need the libsort.so .
  A simple "make" in the package directory should rebuild the 
two libraries from the .c files. 

  If you wish to use the libraries in c code, you should 
update your .bashrc or .cshrc, .whateverrc  to add:
  - the argsort and kdtree directories to CPATH
  - The package directory to LIBRARY_PATH and LD_LIBRARY_PATH

Python code: Julien Dorval
C code: Adapted from c++ to c by Julien Dorval
        original c++ code from Numerical Recipees(3rd edition)

Any questions: dorvaljulien@gmail.com

"""



import numpy as np
import random
import ctypes as C
import matplotlib.pyplot as plt
import time


_lib = np.ctypeslib.load_library('libkdtree.so', 
                                 '/home/observatoire/deparis/Quartz/python/')
double_pointer=C.POINTER(C.c_double)
int_pointer=C.POINTER(C.c_int)


class C_Point(C.Structure):
    pass
C_Point._fields_=[("DIM",C.c_int),("x",double_pointer)]


class Box(C.Structure):
    pass
Box._fields_=[("lo",C_Point),("hi",C_Point),
              ("mom",C.c_int),("dau1",C.c_int),("dau2",C.c_int),
              ("ptlo",C.c_int),("pthi",C.c_int)]

class C_Tree(C.Structure):
    pass
C_Tree._fields_=[("N",C.c_int),("DIM",C.c_int),
                 ("pts",C.POINTER(C_Point)),
                 ("ptindx", int_pointer), ("rptindx",int_pointer),
                 ("nboxes",C.c_int), ("boxes",C.POINTER(Box))  ]

class Tree:
    """

    User manual:
    Create the tree with tree=Tree(x,y,z,..)
    ( works with as many dimensions as you like )
    """

    def __init__(self,*vectors):
        c_tree=KDtree(*vectors)
        self.c_tree=c_tree
        self.N=c_tree.N
        self.DIM=c_tree.DIM
        self.nboxes=c_tree.nboxes
        self.boxes=c_tree.boxes
        self.ptindx=c_tree.ptindx
        self.coords=np.swapaxes([vec for vec in vectors],0,1)

    def draw(self):
        """Visualize a 2d tree"""
        if self.DIM != 2: 
            raise Exception("Cannot draw non 2D tree. Sorry about that.")
        plt.plot(self.coords[:,0],self.coords[:,1],'bo')
        for i in range(self.nboxes): 
            draw_box(self.boxes[i])
        M=1.05*np.max(abs(self.coords))
        plt.xlim(-M,M)
        plt.ylim(-M,M)

    def nearest( self , point ):
        """
        Find the closest neighbour of an arbitrary point x:
            nb = tree.nearest(x)
        """
        cpoint=C_Point()
        cpoint.DIM = len(point)
        cpoint.x= (len(point)*C.c_double)(*point)
        _lib.nearest.argtype = [ C_Tree, C_Point  ]
        _lib.nearest.restype = C.c_int
        return _lib.nearest( self.c_tree, cpoint )

    def nnearest( self, n_point, n_neighbours ):
        """
        Find the n neighbours (identity and distances) of point j:
            indexes, distances = tree.nnearest(j,n)
        """
        _lib.nearest.argtype = [ C_Tree, C.c_int, int_pointer, 
                                 double_pointer, C.c_int ]
        _lib.nnearest.restype = C.c_void_p
        ind_neighbours=(n_neighbours*C.c_int)()
        distances=(n_neighbours*C.c_double)()
        _lib.nnearest(self.c_tree, n_point, 
                      ind_neighbours, distances, n_neighbours)
        ind,dist=[],[]
        for i in range(n_neighbours):
            ind.append(ind_neighbours[i])
            dist.append(distances[i])
        return ind,dist

    def locatenear( self, point, r, nmax=10000000):
        """
        Find all points within r of an arbitrary point x:
             nblist = tree.locatenear(x,r)
        """
        cpoint=C_Point()
        cpoint.DIM = len(point)
        cpoint.x= (len(point)*C.c_double)(*point)
        _lib.locatenear.argtype = [ C_Tree, C_Point, 
                                    C.c_double, int_pointer, C.c_int  ]
        _lib.locatenear.restype = C.c_int
        nb_list=(self.c_tree.N * C.c_int)()
        r = (C.c_double)(r) 
        n_nb=_lib.locatenear(self.c_tree, cpoint, r,nb_list, nmax )
        return np.array(nb_list[:n_nb])





def KDtree( *vectors ):
    """
    This function is the actual kdtree builder, called by the Tree class.
    It accepts as many 1d numpy array as you want dimensions
    """
    _lib.KDtree.argtype = [ double_pointer, C.c_int ]
    _lib.KDtree.restype =  C_Tree
    DIM=len(vectors)
    N=len(vectors[0])
    coord=(DIM*N * C.c_double)()
    for i in range(N):
        for ndim,x in enumerate(vectors):
            coord[N*ndim+i]=x[i]
    return _lib.KDtree(coord,N,DIM)



def draw_box(box):
    xlo,ylo,xhi,yhi=box.lo.x[0],box.lo.x[1],box.hi.x[0],box.hi.x[1],
    plt.plot([xlo,xlo],[ylo,yhi],'k')
    plt.plot([xlo,xhi],[yhi,yhi],'k')
    plt.plot([xhi,xhi],[yhi,ylo],'k')
    plt.plot([xhi,xlo],[ylo,ylo],'k')


"""
if __name__ == "__main__":

    Testing the various functions for a 2d tree.
    
    N=1000000
    
    x = 2*np.random.random(N) -1
    y = 2*np.random.random(N) -1

    fig=plt.figure(1)
    ax=fig.add_subplot(111)
    ax.plot(x,y,'bo')

    t0 = time.time()

    tree = Tree(x,y)

    

    # Testing locatenear:
    if (True):
        print "Testing locatenear"
        point=[0.2,0.2]     # Central point
        radius=0.01         # Maximum distance to the point
        # Rescaling the plot to have a good view of the circle
        zoom=0.5
        plt.xlim((point[0]-(1/zoom)*radius, point[0]+(1/zoom)*radius))
        plt.ylim((point[1]-(1/zoom)*radius, point[1]+(1/zoom)*radius))
        # Plotting the circle around the point:
        th=np.linspace(0,2*np.pi,100)
        xcircle=np.array(point[0]+radius*np.cos(th))
        ycircle=np.array(point[1]+radius*np.sin(th))
        ax.plot(xcircle, ycircle)

        # Plotting central point as a green dot
        ax.plot([point[0],point[0]],[point[1],point[1]],'g.')
        # Getting the indices of all points lying in the circle:
        ind = tree.locatenear(point,radius)
        print ind
        # Plotting these in a different colour:
        for i in ind:
            ax.plot([x[i],x[i]],[y[i],y[i]],'ro')
        plt.show()


    # Testing nnearest:
    if (False):
        print "Testing nnearest"
        index = random.randint(0,N-1)  # getting a particle at random
        Nnb = 4 # How many neighbours we want

        # Getting the indexes and distance to the Nnb closest particles
        ind, dist =tree.nnearest(index, Nnb)
        print "ind :", ind
        print "dist :", dist
        # Rescaling the plot to have a good view of the circle
        zoom=0.5
        maxdist = np.max(dist)
        plt.xlim((x[index]-(1/zoom)*maxdist, x[index]+(1/zoom)*maxdist))
        plt.ylim((y[index]-(1/zoom)*maxdist, y[index]+(1/zoom)*maxdist))
        # Plotting a circle enclosing all the neighbours
        th=np.linspace(0,2*np.pi,100)
        xcircle=np.array(x[index]+maxdist*np.cos(th))
        ycircle=np.array(y[index]+maxdist*np.sin(th))
        ax.plot(xcircle, ycircle)

        # Plotting all neighbours in red
        for i in ind:
            ax.plot([x[i],x[i]],[y[i],y[i]],'ro')
        # Plotting central particle in green
        ax.plot([x[index],x[index]],[y[index],y[index]],'go')
        
        plt.show()
    


    # Testing nearest:
    if (False):
        print "Testing nearest"
        point=[0.2,0.2]
        ind = tree.nearest(point)
        dist= np.sqrt( (point[0]-x[ind])**2 + (point[1]-y[ind])**2  )
        zoom=0.001
        plt.xlim((point[0]-(1/zoom)*dist, point[0]+(1/zoom)*dist))
        plt.ylim((point[1]-(1/zoom)*dist, point[1]+(1/zoom)*dist))
        # Plotting central point as a green dot
        ax.plot([point[0],point[0]],[point[1],point[1]],'g.')
        # Plotting closest particle in red
        ax.plot([x[ind],x[ind]],[y[ind],y[ind]],'ro')
        plt.show()

"""
