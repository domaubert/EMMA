import sys
import os 
import numpy as np
from fonction_IO import *
from fonction_part import *
from fonction_physique import *
	
def PartToVisit(N,t,parts, NameFileOut) :
	print "writing file", NameFileOut

	f = open(NameFileOut , "wb")
	f.write("x y z value\n");


	t=a2t(t)


	inc = 1
	p=0
	while p<N:
		age = t - parts[p].age

		f.write("%g %g %g %g\n" % (parts[p].x,parts[p].y,parts[p].z, age  ))
		p += inc

	f.close()	


def p2v(args, case): 

	filename = args.files[0]

	if   case == 0:
		N,t,part = ReadPart(filename, args)
	elif case == 1:
		N,t,part = ReadStars(filename, args)
		filename = filename[:-17] +  "star" +  filename[-13:]

	if N :	
		PartToVisit(N,t,part, filename + ".3D")
	else :
		print "pas de particules"



def FolderToVisit(foldername):

	files = listPart(foldername)

	for i in range(len(files)):
		file = files[i]
		if file[0:3]== "par" and file[-1]!="D":	
			N,t,parts= ReadPart(file)		
			PartToVisit(parts, file+".3D")
#			PartToVisit(getStars(parts), foldername + file +".stars.3D")


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


