import sys, os
import argparse
import numpy as np
from struct import *

class Param : 
	def __init__(self,filename):

		file = open(filename, "rb")

		self.PIC          = unpack('i',file.read(4))[0]
		self.WHYDRO2      = unpack('i',file.read(4))[0]
		self.WGRAV        = unpack('i',file.read(4))[0]
		self.WRAD         = unpack('i',file.read(4))[0]
		self.WRADHYD      = unpack('i',file.read(4))[0]
		self.TESTCOSMO    = unpack('i',file.read(4))[0]
		self.WDBG         = unpack('i',file.read(4))[0]
		self.STARS        = unpack('i',file.read(4))[0]
		self.WMPI         = unpack('i',file.read(4))[0]
		self.FLOORDT      = unpack('i',file.read(4))[0]
		self.WCUDA_ERR    = unpack('i',file.read(4))[0]
		self.NOCOMP       = unpack('i',file.read(4))[0]
		self.GRAFIC       = unpack('i',file.read(4))[0]
		self.ZELDOVICH    = unpack('i',file.read(4))[0]
		self.EVRARD       = unpack('i',file.read(4))[0]
		self.EDBERT       = unpack('i',file.read(4))[0]
		self.TUBE         = unpack('i',file.read(4))[0]
		self.PARTN        = unpack('i',file.read(4))[0]
		self.PART2        = unpack('i',file.read(4))[0]
		self.WRADTEST     = unpack('i',file.read(4))[0]
		self.TESTCLUMP    = unpack('i',file.read(4))[0]
		self.PART_EGY     = unpack('i',file.read(4))[0]
		self.PERFECT      = unpack('i',file.read(4))[0]
		self.FASTGRAV     = unpack('i',file.read(4))[0]
		self.ONFLYRED     = unpack('i',file.read(4))[0]
		self.RIEMANN_HLLC = unpack('i',file.read(4))[0]
		self.RIEMANN_EXACT= unpack('i',file.read(4))[0]
		self.PRIMITIVE    = unpack('i',file.read(4))[0]
		self.DUAL_E       = unpack('i',file.read(4))[0]
		self.WCHEM        = unpack('i',file.read(4))[0]
		self.S_100000     = unpack('i',file.read(4))[0]
		self.COOLING      = unpack('i',file.read(4))[0]
		self.UVBKG        = unpack('i',file.read(4))[0]
		self.TRANSZM      = unpack('i',file.read(4))[0]
		self.TRANSZP      = unpack('i',file.read(4))[0]
		self.TRANSYM      = unpack('i',file.read(4))[0]
		self.TRANSYP      = unpack('i',file.read(4))[0]
		self.TRANSXM      = unpack('i',file.read(4))[0]
		self.TRANSXP      = unpack('i',file.read(4))[0]
		self.REFXM        = unpack('i',file.read(4))[0]
		self.REFYM        = unpack('i',file.read(4))[0]
		self.REFZM	  = unpack('i',file.read(4))[0]

		self.npartmax     = unpack('i',file.read(4))[0]
		self.ngridmax     = unpack('i',file.read(4))[0]
		self.nbuff        = unpack('i',file.read(4))[0]
		self.ndumps       = unpack('i',file.read(4))[0]
		self.nsteps       = unpack('i',file.read(4))[0]

		self.lcoarse      = unpack('i',file.read(4))[0]
		self.lmax         = unpack('i',file.read(4))[0]

		self.niter        = unpack('i',file.read(4))[0]

		self.gstride      = unpack('i',file.read(4))[0]
		self.hstride      = unpack('i',file.read(4))[0]

		self.dt           = unpack('f',file.read(4))[0]
		self.tmax         = unpack('f',file.read(4))[0]
		self.time_max     = unpack('f',file.read(4))[0]

		self.maxhash      = unpack('i',file.read(4))[0]

		self.amrthresh    = unpack('f',file.read(4))[0]
		self.nsmooth      = unpack('i',file.read(4))[0]

		self.poissonacc   = unpack('f',file.read(4))[0]
		self.mgridlmin    = unpack('i',file.read(4))[0]
		self.nvcycles     = unpack('i',file.read(4))[0]
		self.nrelax       = unpack('i',file.read(4))[0]

		self.nrestart     = unpack('i',file.read(4))[0]
		self.nsubcycles   = unpack('i',file.read(4))[0]

		self.nthread      = unpack('i',file.read(4))[0]
		self.nstream      = unpack('i',file.read(4))[0]
		
		self.egy_rhs      = unpack('f',file.read(4))[0]
		self.egy_0        = unpack('f',file.read(4))[0]
		self.egy_last     = unpack('f',file.read(4))[0]
		self.egy_timelast = unpack('f',file.read(4))[0]
		self.egy_totlast  = unpack('f',file.read(4))[0]
		
		if self.WRAD :
			self.unit_l                = unpack('f',file.read(4))[0]
			self.unit_v                = unpack('f',file.read(4))[0]
			self.unit_t                = unpack('f',file.read(4))[0]
			self.unit_n                = unpack('f',file.read(4))[0]
			self.unit_mass             = unpack('f',file.read(4))[0]

			self.clight                = unpack('f',file.read(4))[0] 	
			self.fudgecool             = unpack('f',file.read(4))[0] 	
			self.ncvgcool              = unpack('i',file.read(4))[0] 	
	  
			self.denthresh             = unpack('f',file.read(4))[0]	
			self.tmpthresh             = unpack('f',file.read(4))[0]	
			self.srcint                = unpack('f',file.read(4))[0]	

		if self.TESTCOSMO :
			self.om                    = unpack('f',file.read(4))[0]
			self.ov                    = unpack('f',file.read(4))[0]
			self.ob                    = unpack('f',file.read(4))[0]
			self.H0                    = unpack('f',file.read(4))[0]

		if self.STARS :
			self.overdensity_cond      = unpack('f',file.read(4))[0]
			self.density_cond          = unpack('f',file.read(4))[0]
			self.tcar                  = unpack('f',file.read(4))[0]
		

		file.close()

def readFile(filename):

	try:
		f = open(filename, "r")	
	except IOError:
		print 'cannot open', filename 

	data = f.read()
	f.close()
	return data	


def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False


def getSnapMax(args):
	
	print args.folder[0]
	files = sorted(os.listdir(args.folder[0]))

	idx = len(files)-1
	while 1 :
		f = files[idx]
		if f[0:4]== "part" or f[0:4]== "grid" or f[0:4]== "star":	
			break
		idx -= 1

	args.nmax = int(f[-12:-7])

def listPart(foldername):

	files = os.listdir(foldername)

	tmp=[]
	for file in files:
		if file[0:4]== "part" :	
			tmp.append(foldername + file)

	return  sorted(tmp)


def listOriginePart(foldername) :

	files = listPart(foldername)

	tmp=[]
	for file in files :
		if file[-3:]!= ".3D" and file[-7:]== ".p00000" :
			tmp.append( file[:-7])	
	return tmp


def getNsnap(foldername):
	return len(listOriginePart(foldername))-1

def getNproc(args):
	files = os.listdir(args.folder[0])
	tmp =0
	for file in files:
		if file[0:10]=="part.00000" and file[-3:]!=".3D" and file[-5:]!=".cube":	
			tmp +=1
	args.nproc = tmp


def snaprange(args):
	
	if args.fro != -1 or args.to != -1 : 
		if args.fro != -1 : 
			fro = args.fro
		else :
			fro = 0

		if args.to != -1 : 
			to = args.to
		else :
			to = getNsnap(args.folder[0])

		args.snap = np.arange(fro, to+1)

def num2snap(args):
	args.files = listOriginePart(args.folder[0])
	if args.snap != []:
		f=[]
		for i in args.snap:
			f.append(args.files[i])
		args.files = f

def getA(filename):
	filePart = open(filename, "rb")
	filePart.read(4)
	a = unpack('f',filePart.read(4))[0]
	filePart.close()
	return a

def getZ(filename):
	return 1.0/getA(filename) -1

def getFolder(args):
	for i in range(len(args.files)):
		args.folder.append(args.files[0][:-10])

def getargs() :
	parser = argparse.ArgumentParser(description='EMMA analyse program')

	parser.add_argument('-fo',   action="append",     default=[],            help = "witch folder to use", dest = "folder")
	parser.add_argument('-fi',   action="append",     default=[],            help = "snapshot file(s)", dest = "files")
	parser.add_argument('-n',    action="append",     default=[], type=int,  help = "snapshot number n", dest = "snap")
	parser.add_argument('-from', action="store",      default=-1, type=int,  help = "snapshot number to begin with, 0 by default", dest = "fro")
	parser.add_argument('-to',   action="store",      default=-1, type=int,  help = "snapshot number to end with, end by default", dest = "to")
	parser.add_argument('-plot', action="store_true", default= False,        help = "plot the current set of values")
	parser.add_argument('-part', action="store_true", default= True,         help = "is it a particle snap?")
	parser.add_argument('-grid', action="store_true", default= False,        help = "is it a grid snap?")
	parser.add_argument('-np',   action="store",      default=-1, type=int,  help = "number of procs used to generate the snap. only usefull to force it", dest = "nproc")
	parser.add_argument('-l'    ,action="store",      default=0 , type=int,  help = "level param of oct2grid", dest = "level")
	parser.add_argument('-field',action="append",     default=["field.d"] ,  help = "field param of oct2grid", dest = "field")

	parser.add_argument('-tag',  action="append",     default=[], type=int,  help = "halo number", dest = "tag")

	parser.add_argument('-nmax', action="store",      default=-1, type=int,  help = "max snapshot number", dest = "nmax")
	parser.add_argument('-nmin', action="store",      default=0, type=int,   help = "min snapshot number", dest = "nmin")



	args = parser.parse_args()




	if args.folder == []:
		if args.files !=[] :
			getFolder(args)
		else :
			args.folder.append("data/")

	for fol in args.folder:
		if (fol[-1]!="/"):
			fol += "/"




	if args.nmax == -1:
		getSnapMax(args)

	if args.nproc == -1:
		nproc = getNproc(args)





	snaprange(args)

	if args.files == []:
		num2snap(args)
	

	return args







