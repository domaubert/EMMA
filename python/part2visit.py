#!/usr/bin/env python
from fonction_visit import *
from fonction_part import *


if __name__ == "__main__":
	


	args=getargs()
	foldername=args.folder[0]

	if args.files != []:
		filename = args.files[0]

		'''
		N,t,parts= ReadPart(filename, args)

		if N :	
			PartToVisit(parts, filename + ".3D")
		else :
			print "pas de particules"
		'''

		Nstars,t,stars = ReadStars(filename, args)

		if Nstars :	
			filename = filename[:-17] +  "star" +  filename[-13:]
			PartToVisit(Nstars,stars, filename + ".3D")
		else :
			print "pas de particules"


	else:
		FolderToVisit(foldername)
		makevideofile(foldername)





