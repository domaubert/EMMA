import sys
import os 
import numpy as np
from fonction_IO import *

def makevideofile(foldername) :
	
	files = listPart(foldername)

	f = open("../data/stars.visit", "wt")

	tmp=[]
	for file in files:
		if file[0:3]== "par" :	
			tmp.append(file)

	for i in range(len(files)):
		if tmp[i][-9:]==".stars.3D" :
			f.write(files[i] )
			f.write( "\n" )
	f.close()	

	
def PartToVisit(parts, NameFileOut) :
	print "writing file", NameFileOut
	f = open(NameFileOut , "wt")
	f.write("x y z value\n");
	for p in parts:
		if  p.mass < pow(10,-10.25):
			f.write("%g %g %g %g\n" % (p.x,p.y,p.z,p.mass))
	f.close()	



def FolderToVisit(foldername):

	files = listPart(foldername)

	for i in range(len(files)):
		file = files[i]
		if file[0:3]== "par" and file[-1]!="D":	
			N,t,parts= ReadPart(file)		
			PartToVisit(parts, file+".3D")
#			PartToVisit(getStars(parts), foldername + file +".stars.3D")

