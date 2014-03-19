#!/usr/bin/env python
from fonction_visit import *
from fonction_part import *


if __name__ == "__main__":
	


	args=getargs()
	foldername=args.folder[0]

	if args.files != []:

		filename = args.files[0]

		N,t,parts= ReadPart(filename, args)
		Nstars, stars = getStars(parts)

		PartToVisit(parts, filename + ".3D")	
		PartToVisit(stars, filename +".stars.3D")
		
	else:
		FolderToVisit(foldername)
		makevideofile(foldername)





