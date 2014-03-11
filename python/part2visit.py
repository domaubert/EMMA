from fonction_visit import *
from fonction_part import *


if __name__ == "__main__":
	


	args=getargs()
	foldername=args.folder

	if args.files != []:

		filename = args.files[0]

		N,t,parts= ReadPart(filename, args)
		Nstars, stars = getStars(parts)

		PartToVisit(stars, filename +".stars.3D")



		'''
		if sys.argv[1]=="parts" :
			PartToVisit(parts, filename + ".3D")

		elif sys.argv[1]=="stars" :
			PartToVisit(getStars(parts), filename +".stars.3D")

		elif sys.argv[1]=="all" :
			PartToVisit(parts, filename + ".3D")
			PartToVisit(getStars(parts), filename +".stars.3D")
		else :
			print "error"
		'''
	else:
		FolderToVisit(foldername)
		makevideofile(foldername)





