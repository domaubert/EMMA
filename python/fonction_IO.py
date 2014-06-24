import sys, os
import argparse
import numpy as np
from struct import *


def inc(x):
	x += 1 
	return x 


class Param : 
	def __init__(self,folder):

		file = open(folder + "param.00000.p00000", "rb")
		data = file.read()
		data = [x for x in data.split('\n')]
		for i in range(len(data)):
			data[i] = [x for x in data[i].split('\t')] 

		i =0

		self.PIC          = int(data[i][1]); i += 1
		self.WHYDRO2      = int(data[i][1]); i += 1
		self.WGRAV        = int(data[i][1]); i += 1
		self.WRAD         = int(data[i][1]); i += 1
		self.WRADHYD      = int(data[i][1]); i += 1
		self.TESTCOSMO    = int(data[i][1]); i += 1
		self.WDBG         = int(data[i][1]); i += 1
		self.STARS        = int(data[i][1]); i += 1
		self.WMPI         = int(data[i][1]); i += 1
		self.FLOORDT      = int(data[i][1]); i += 1
		self.WCUDA_ERR    = int(data[i][1]); i += 1
		self.NOCOMP       = int(data[i][1]); i += 1
		self.GRAFIC       = int(data[i][1]); i += 1
		self.ZELDOVICH    = int(data[i][1]); i += 1
		self.EVRARD       = int(data[i][1]); i += 1
		self.EDBERT       = int(data[i][1]); i += 1
		self.TUBE         = int(data[i][1]); i += 1
		self.PARTN        = int(data[i][1]); i += 1
		self.PART2        = int(data[i][1]); i += 1
		self.WRADTEST     = int(data[i][1]); i += 1
		self.TESTCLUMP    = int(data[i][1]); i += 1
		self.PART_EGY     = int(data[i][1]); i += 1
		self.PERFECT      = int(data[i][1]); i += 1
		self.FASTGRAV     = int(data[i][1]); i += 1
		self.ONFLYRED     = int(data[i][1]); i += 1
		self.RIEMANN_HLLC = int(data[i][1]); i += 1
		self.RIEMANN_EXACT= int(data[i][1]); i += 1
		self.PRIMITIVE    = int(data[i][1]); i += 1
		self.DUAL_E       = int(data[i][1]); i += 1
		self.WCHEM        = int(data[i][1]); i += 1
		self.S_100000     = int(data[i][1]); i += 1
		self.COOLING      = int(data[i][1]); i += 1
		self.UVBKG        = int(data[i][1]); i += 1
		self.TRANSZM      = int(data[i][1]); i += 1
		self.TRANSZP      = int(data[i][1]); i += 1
		self.TRANSYM      = int(data[i][1]); i += 1
		self.TRANSYP      = int(data[i][1]); i += 1
		self.TRANSXM      = int(data[i][1]); i += 1
		self.TRANSXP      = int(data[i][1]); i += 1
		self.REFXM        = int(data[i][1]); i += 1
		self.REFYM        = int(data[i][1]); i += 1
		self.REFZM	  = int(data[i][1]); i += 1

		self.npartmax     = int(data[i][1]); i += 1
		self.ngridmax     = int(data[i][1]); i += 1
		self.nbuff        = int(data[i][1]); i += 1
		self.ndumps       = int(data[i][1]); i += 1
		self.nsteps       = int(data[i][1]); i += 1

		self.lcoarse      = int(data[i][1]); i += 1
		self.lmax         = int(data[i][1]); i += 1

		self.niter        = int(data[i][1]); i += 1

		self.gstride      = int(data[i][1]); i += 1
		self.hstride      = int(data[i][1]); i += 1

		self.dt           = float(data[i][1]); i += 1
		self.tmax         = float(data[i][1]); i += 1
		self.time_max     = float(data[i][1]); i += 1

		self.maxhash      = int(data[i][1]); i += 1

		self.amrthresh    = float(data[i][1]); i += 1
		self.nsmooth      = int(data[i][1]); i += 1

		self.poissonacc   = float(data[i][1]); i += 1
		self.mgridlmin    = int(data[i][1]); i += 1
		self.nvcycles     = int(data[i][1]); i += 1
		self.nrelax       = int(data[i][1]); i += 1

		self.nrestart     = int(data[i][1]); i += 1
		self.nsubcycles   = int(data[i][1]); i += 1

		self.nthread      = int(data[i][1]); i += 1
		self.nstream      = int(data[i][1]); i += 1
		
		self.egy_rhs      = float(data[i][1]); i += 1
		self.egy_0        = float(data[i][1]); i += 1
		self.egy_last     = float(data[i][1]); i += 1
		self.egy_timelast = float(data[i][1]); i += 1
		self.egy_totlast  = float(data[i][1]); i += 1
		
		if self.WRAD :
			self.unit_l                = float(data[i][1]); i += 1
			self.unit_v                = float(data[i][1]); i += 1
			self.unit_t                = float(data[i][1]); i += 1
			self.unit_n                = float(data[i][1]); i += 1
			self.unit_mass             = float(data[i][1]); i += 1

			self.clight                = float(data[i][1]); i += 1 	
			self.fudgecool             = float(data[i][1]); i += 1 	
			self.ncvgcool              = int(data[i][1]); i += 1 	
	  
			self.denthresh             = float(data[i][1]); i += 1	
			self.tmpthresh             = float(data[i][1]); i += 1	
			self.srcint                = float(data[i][1]); i += 1	

		if self.TESTCOSMO :
			self.om                    = float(data[i][1]); i += 1
			self.ov                    = float(data[i][1]); i += 1
			self.ob                    = float(data[i][1]); i += 1
			self.H0                    = float(data[i][1]); i += 1

		if self.STARS :
			self.overdensity_cond      = float(data[i][1]); i += 1
			self.density_cond          = float(data[i][1]); i += 1
			self.tcar                  = float(data[i][1]); i += 1
			self.tlife    		   = float(data[i][1]); i += 1
			self.feedback_eff          = float(data[i][1]); i += 1
			self.feedback_frac         = float(data[i][1]); i += 1
		

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


def getSlurm(args):

	files = os.listdir(args.folder[0])

	for file in files:
		if file[-19:-14]== "slurm" : 	
			return args.folder[0] + file




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

def getNp(name):

	files = os.listdir(name[:-10])

	tmp =0
	for file in files:
		if name[-10:] in file :
			tmp +=1
	return tmp


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
			print args.folder
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
	


	P = Param(args.folder[0])

	return args







