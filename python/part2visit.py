
from fonction_part import *



if __name__ == "__main__":
	

	foldername="../data/"

	if len(sys.argv)>2 :

		filename = foldername + "part." + sys.argv[2].zfill(5) + ".p00000"
		N,t,parts= ReadPart(filename)

		if sys.argv[1]=="parts" :
			PartToVisit(parts, filename + ".3D")

		elif sys.argv[1]=="stars" :
			PartToVisit(getStars(parts), filename +".stars.3D")

		elif sys.argv[1]=="all" :
			PartToVisit(parts, filename + ".3D")
			PartToVisit(getStars(parts), filename +".stars.3D")
		else :
			print "error"

	else:
		FolderToVisit(foldername)
		makevideofile(foldername)





