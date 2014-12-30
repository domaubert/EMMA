#!/usr/bin/env python
from fonction_visit import *
from fonction_part import*

if __name__ == "__main__":
	


	args=getargs()


	if args.files != []:

		case = 1	# 0->Part	1->Stars	2->both

		if   case == 0 : 
			p2v(args, 0)	
		elif case == 1:
			p2v(args, 1)	
		elif case == 2:
			p2v(args, 0)
			p2v(args, 1)


	else:
		FolderToVisit(foldername)
		makevideofile(foldername)





