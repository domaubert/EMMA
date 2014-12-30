import sys,os
import numpy as np
import argparse

class Params:
	def __init__(self,folder = "data"):
		filename = "%sparam.00000.p00000"%folder
		self.d = {}
		with open(filename) as f:
			for line in f:
				(key, val) = line.split()
				self.d[key] = val					

	def get(self):
		return self.d

def getNproc(filepath="data/grid.00000"):
	folder, filename = filepath.split("/")
	files = os.listdir(folder)
	nProc =0
	for file in files:
		if filename in file :
			nProc += 1
	return nProc
	
def getargs() :
	parser = argparse.ArgumentParser(description='EMMA analyse program')

	parser.add_argument('-fo',   action="store",      default=["data/"],     help = "witch folder to use", dest = "folder")
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
	return args
