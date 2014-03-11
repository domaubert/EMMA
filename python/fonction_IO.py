import sys
import os 
import argparse
import numpy as np
import matplotlib.pylab as plt

class Param : 
	def __init__(self,file):


		self.PIC          = np.fromfile(file, dtype=np.int32,   count=1)
		self.WHYDRO2      = np.fromfile(file, dtype=np.int32,   count=1)
		self.WGRAV        = np.fromfile(file, dtype=np.int32,   count=1)
		self.WRAD         = np.fromfile(file, dtype=np.int32,   count=1)
		self.WRADHYD      = np.fromfile(file, dtype=np.int32,   count=1)
		self.TESTCOSMO    = np.fromfile(file, dtype=np.int32,   count=1)
		self.WDBG         = np.fromfile(file, dtype=np.int32,   count=1)
		self.STARS        = np.fromfile(file, dtype=np.int32,   count=1)
		self.WMPI         = np.fromfile(file, dtype=np.int32,   count=1)
		self.FLOORDT      = np.fromfile(file, dtype=np.int32,   count=1)
		self.WCUDA_ERR    = np.fromfile(file, dtype=np.int32,   count=1)
		self.NOCOMP       = np.fromfile(file, dtype=np.int32,   count=1)
		self.GRAFIC       = np.fromfile(file, dtype=np.int32,   count=1)
		self.ZELDOVICH    = np.fromfile(file, dtype=np.int32,   count=1)
		self.EVRARD       = np.fromfile(file, dtype=np.int32,   count=1)
		self.EDBERT       = np.fromfile(file, dtype=np.int32,   count=1)
		self.TUBE         = np.fromfile(file, dtype=np.int32,   count=1)
		self.PARTN        = np.fromfile(file, dtype=np.int32,   count=1)
		self.PART2        = np.fromfile(file, dtype=np.int32,   count=1)
		self.WRADTEST     = np.fromfile(file, dtype=np.int32,   count=1)
		self.TESTCLUMP    = np.fromfile(file, dtype=np.int32,   count=1)
		self.PART_EGY     = np.fromfile(file, dtype=np.int32,   count=1)
		self.PERFECT      = np.fromfile(file, dtype=np.int32,   count=1)
		self.FASTGRAV     = np.fromfile(file, dtype=np.int32,   count=1)
		self.ONFLYRED     = np.fromfile(file, dtype=np.int32,   count=1)
		self.RIEMANN_HLLC = np.fromfile(file, dtype=np.int32,   count=1)
		self.RIEMANN_EXACT= np.fromfile(file, dtype=np.int32,   count=1)
		self.PRIMITIVE    = np.fromfile(file, dtype=np.int32,   count=1)
		self.DUAL_E       = np.fromfile(file, dtype=np.int32,   count=1)
		self.WCHEM        = np.fromfile(file, dtype=np.int32,   count=1)
		self.S_100000     = np.fromfile(file, dtype=np.int32,   count=1)
		self.COOLING      = np.fromfile(file, dtype=np.int32,   count=1)
		self.UVBKG        = np.fromfile(file, dtype=np.int32,   count=1)
		self.TRANSZM      = np.fromfile(file, dtype=np.int32,   count=1)
		self.TRANSZP      = np.fromfile(file, dtype=np.int32,   count=1)
		self.TRANSYM      = np.fromfile(file, dtype=np.int32,   count=1)
		self.TRANSYP      = np.fromfile(file, dtype=np.int32,   count=1)
		self.TRANSXM      = np.fromfile(file, dtype=np.int32,   count=1)
		self.TRANSXP      = np.fromfile(file, dtype=np.int32,   count=1)
		self.REFXM        = np.fromfile(file, dtype=np.int32,   count=1)
		self.REFYM        = np.fromfile(file, dtype=np.int32,   count=1)
		self.REFZM	  = np.fromfile(file, dtype=np.int32,   count=1)

		self.npartmax     = np.fromfile(file, dtype=np.int32,   count=1)
		self.ngridmax     = np.fromfile(file, dtype=np.int32,   count=1)
		self.nbuff        = np.fromfile(file, dtype=np.int32,   count=1)
		self.ndumps       = np.fromfile(file, dtype=np.int32,   count=1)
		self.nsteps       = np.fromfile(file, dtype=np.int32,   count=1)

		self.lcoarse      = np.fromfile(file, dtype=np.int32,   count=1)
		self.lmax         = np.fromfile(file, dtype=np.int32,   count=1)

		self.niter        = np.fromfile(file, dtype=np.int32,   count=1)

		self.gstride      = np.fromfile(file, dtype=np.int32,   count=1)
		self.hstride      = np.fromfile(file, dtype=np.int32,   count=1)

		self.dt           = np.fromfile(file, dtype=np.float32, count=1)
		self.tmax         = np.fromfile(file, dtype=np.float32, count=1)
		self.time_max     = np.fromfile(file, dtype=np.float32, count=1)

		self.maxhash      = np.fromfile(file, dtype=np.int32,   count=1)

		self.amrthresh    = np.fromfile(file, dtype=np.float32, count=1)
		self.nsmooth      = np.fromfile(file, dtype=np.int32,   count=1)

		self.poissonacc   = np.fromfile(file, dtype=np.float32, count=1)
		self.mgridlmin    = np.fromfile(file, dtype=np.int32,   count=1)
		self.nvcycles     = np.fromfile(file, dtype=np.int32,   count=1)
		self.nrelax       = np.fromfile(file, dtype=np.int32,   count=1)

		self.nrestart     = np.fromfile(file, dtype=np.int32,   count=1)
		self.nsubcycles   = np.fromfile(file, dtype=np.int32,   count=1)

		self.nthread      = np.fromfile(file, dtype=np.int32,   count=1)
		self.nstream      = np.fromfile(file, dtype=np.int32,   count=1)
		
		self.egy_rhs      = np.fromfile(file, dtype=np.float32, count=1)
		self.egy_0        = np.fromfile(file, dtype=np.float32, count=1)
		self.egy_last     = np.fromfile(file, dtype=np.float32, count=1)
		self.egy_timelast = np.fromfile(file, dtype=np.float32, count=1)
		self.egy_totlast  = np.fromfile(file, dtype=np.float32, count=1)
		
		if self.WRAD :
			self.unit_l                = np.fromfile(file, dtype=np.float32, count=1)
			self.unit_v                = np.fromfile(file, dtype=np.float32, count=1)
			self.unit_t                = np.fromfile(file, dtype=np.float32, count=1)
			self.unit_n                = np.fromfile(file, dtype=np.float32, count=1)
			self.unit_mass             = np.fromfile(file, dtype=np.float32, count=1)

			self.clight                = np.fromfile(file, dtype=np.float32, count=1) 	
			self.fudgecool             = np.fromfile(file, dtype=np.float32, count=1) 	
			self.ncvgcool              = np.fromfile(file, dtype=np.int32,   count=1) 	
	  
			self.denthresh             = np.fromfile(file, dtype=np.float32, count=1)	
			self.tmpthresh             = np.fromfile(file, dtype=np.float32, count=1)	
			self.srcint                = np.fromfile(file, dtype=np.float32, count=1)	

		if self.TESTCOSMO :
			self.om                    = np.fromfile(file, dtype=np.float32, count=1)
			self.ov                    = np.fromfile(file, dtype=np.float32, count=1)
			self.ob                    = np.fromfile(file, dtype=np.float32, count=1)
			self.H0                    = np.fromfile(file, dtype=np.float32, count=1)

		if self.STARS :
			self.overdensity_cond      = np.fromfile(file, dtype=np.float32, count=1)
			self.density_cond          = np.fromfile(file, dtype=np.float32, count=1)
			self.t_car                 = np.fromfile(file, dtype=np.float32, count=1)
			self.eff                   = np.fromfile(file, dtype=np.float32, count=1)
	


def listPart(foldername):

	files = os.listdir(foldername)

	tmp=[]
	for file in files:
		if file[0:3]== "par" :	
			tmp.append(foldername + file)

	return  sorted(tmp)

def listOriginePart(foldername) :

	files = listPart(foldername)

	tmp=[]
	for file in files :
		if file[-3:]!= ".3D" and file[-6]== "p"  and file[-1]== "0" :
			tmp.append( file)
			

	return tmp

def getNsnap(foldername):
	return len(listOriginePart(foldername))-1

def getNproc(args):
	files = os.listdir(args.folder)
	tmp =0
	for file in files:
		if file[0:10]=="part.00000" and file[-3:]!=".3D" :	
			tmp +=1
	return tmp


def snaprange(args):
	
	if args.fro != -1 and args.to != -1 : 
		if args.fro != -1 : 
			fro = args.fro
		else :
			fro = 0

		if args.to != -1 : 
			to = args.to
		else :
			to = getNsnap(args.folder)

		args.snap = np.arange(fro, to+1)

def num2snap(args):
	args.files = listOriginePart(args.folder)
	if args.snap != []:
		f=[]
		for i in args.snap:
			f.append(args.files[i])
		args.files = f


def getargs() :
	parser = argparse.ArgumentParser(description='EMMA analyse program')

	parser.add_argument('-fo',   action="store",      default="../data/",    help = "witch folder to use", dest = "folder")
	parser.add_argument('-fi',   action="append",     default=[],            help = "snapshot file(s)", dest = "files")
	parser.add_argument('-n',    action="append",     default=[], type=int,  help = "snapshot number n", dest = "snap")
	parser.add_argument('-from', action="store",      default=-1, type=int,  help = "snapshot number to begin with, 0 by default", dest = "fro")
	parser.add_argument('-to',   action="store",      default=-1, type=int,  help = "snapshot number to end with, end by default", dest = "to")
	parser.add_argument('-plot', action="store_true", default= False,        help = "plot the current set of values")
	parser.add_argument('-part', action="store_true", default= True,         help = "is it a particle snap?")
	parser.add_argument('-grid', action="store_true", default= False,        help = "is it a grid snap?")
	parser.add_argument('-np',   action="store",      default=[], type=int,  help = "number of procs used to generate the snap. only usefull to force it", dest = "nproc")

	args = parser.parse_args()

	
	snaprange(args)
	num2snap(args)




	return args
	

if __name__ == "__main__":
	'''
	file = open("../data/param", "rb")
	p = Param(file)
	print p.nstream
	'''
	args = getargs()
	print(args.snap)
	print(args.to)
