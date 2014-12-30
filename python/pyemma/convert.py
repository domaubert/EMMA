import sys, os
import numpy as np

from struct import *

from part import *
from physique import *
from amr import *

def star2obj(filename):
	"""
		Convert EMMA star file to one standard .obj file.
	"""

	N,a,parts = readPart(filename)
	t = a2t(a)

	scale = 64

	file = open("star.obj", "w")

	for i in range(N):
		x = (parts.x[i]-0.5)*scale
		y = (parts.y[i]-0.5)*scale
		z = (parts.z[i]-0.5)*scale

		file.write("v\t%f\t%f\t%f\n"%( x,y,z ))

	file.flush()
	file.close()

def sortStar(filename):
	"""
		Separate the three type of star :  	active, supernovae and dead 
		fonction of there age.
		active if age < 3e7 yr
		supernovae if 3e7 < age < 3e8 yr
		dead if 3e8 < age
	"""
	N,a,parts = readPart(filename)
	t = a2t(a)

	scale = 64

	active_star = []
	supernovae  = []
	dead_star = []

	for i in range(N):
		age = t - parts.age[i]
		x = (parts.x[i]-0.5)*scale
		y = (parts.y[i]-0.5)*scale
		z = (parts.z[i]-0.5)*scale

		if age < 3e7: 
			active_star.append((x,y,z))
		elif age < 3e8: 
			supernovae.append((x,y,z))
		else : 
			dead_star.append((x,y,z))

	return active_star, supernovae, dead_star


def star2obj_v2(filename):
	"""
		Convert EMMA star file to one standard .obj 
		and separate the three type of star :  	active, supernovae and dead 
		to three sub object

		!!!!!!!!!!!!!!!!!!!!!!!!
		(dont work right now)
		!!!!!!!!!!!!!!!!!!!!!!!!
	"""
 	active_star, supernovae, dead_star =	sortStar()

	file = open("star.obj", "w")

	file.write("# cube.obj\n\n")

	file.write("o active\n\n")
	for i in range(len(active_star)):
			file.write("v\t%f\t%f\t%f\n"%( active_star[i][0],active_star[i][1],active_star[i][2] ))

	file.write("\no supernovae\n\n")
	for i in range(len(supernovae)):
			file.write("v\t%f\t%f\t%f\n"%( supernovae[i][0],supernovae[i][1],supernovae[i][2] ))

	file.write("\no dead\n\n")
	for i in range(len(dead_star)):
			file.write("v\t%f\t%f\t%f\n"%( dead_star[i][0],dead_star[i][1],dead_star[i][2] ))

	file.flush()
	file.close()


def star23obj(filename):
	"""
		Convert EMMA star file to three standard .obj 
		by separate the three type of star : active, supernovae and dead 
		to three file
	"""

 	active_star, supernovae, dead_star =	sortStar()

	file_active_star = open("active.obj", "w")
	file_supernovae = open("supernovae.obj", "w")
	file_dead_star = open("dead.obj", "w")

	for i in range(len(active_star)):
			file_active_star.write("v\t%f\t%f\t%f\n"%( active_star[i][0],active_star[i][1],active_star[i][2] ))
	for i in range(len(supernovae)):
			file_supernovae.write("v\t%f\t%f\t%f\n"%( supernovae[i][0],supernovae[i][1],supernovae[i][2] ))
	for i in range(len(dead_star)):
			file_dead_star.write("v\t%f\t%f\t%f\n"%( dead_star[i][0],dead_star[i][1],dead_star[i][2] ))

	file_active_star.flush()
	file_supernovae.flush()
	file_dead_star.flush()

	file_active_star.close()
	file_supernovae.close()
	file_dead_star.close()

def cube2blenderVoxel(filename):
	field = "field.d"
	level = 7

	oct2grid(filename, field, level)
	data = array(filename.replace("grid", "cube"))

	size = data.getSize()
	bin = data.getData()
	header = np.array([size[0],size[1],size[2],1])


	file = open("cube.bvox", "wb")
	header.astype('<i4').tofile(file)
	bin.astype('<f4').tofile(file)
	file.close()

def field2ID(field="field.d"):
	"""
	Convert field from emma description to oct2grid description
	"""

	if field=="density" : 
		f=" -1 "
	if field=="level" : 
		f=" 0 "
	if field=="den" : 
		f=" 1 "
	if field=="pot" : 
		f=" 2 "
	if field=="cpu" : 
		f=" 3 "
	if field=="marked" : 
		f=" 4 "
	if field=="res" : 
		f=" 6 "

	if field=="field.d" : 
		f=" 101 "
	if field=="field.u" : 
		f=" 102 "
	if field=="field.v" : 
		f=" 103 "
	if field=="field.w" : 
		f=" 104 "
	if field=="field.p" : 
		f=" 105 "
	if field=="field.E" : 
		f=" 106 "
		
	if field=="f0" : 
		f=" 201 "
	if field=="f1" : 
		f=" 202 "
	if field=="f2" : 
		f=" 203 "

	if field=="e0" : 
		f=" 701 "
	if field=="fx0" : 
		f=" 702 "		
	if field=="fy0" : 
		f=" 703 "		
	if field=="fz0" : 
		f=" 704 "				
	if field=="e1" : 
		f=" 711 "
	if field=="fx1" : 
		f=" 712 "		
	if field=="fy1" : 
		f=" 713 "		
	if field=="fz1" : 
		f=" 714 "		
	if field=="e2" : 
		f=" 721 "
	if field=="fx2" : 
		f=" 722 "		
	if field=="fy2" : 
		f=" 723 "		
	if field=="fz2" : 
		f=" 724 "	

	if field=="rfield.src" : 
		f=" 705 "
	if field=="rfield.xion" : 
		f=" 706 "
	if field=="rfield.temp" : 
		f=" 707 "
	if field=="snfb" : 
		f=" 709 "
	if field=="eint" : 
		f=" 721 "
	if field=="dX" : 
		f=" 1006 "
	
	return f


def oct2grid(data,level, force=0, nproc=getNproc(), field="field.d", xmin=0, xmax=-1, ymin=0, ymax=-1, zmin=0, zmax=-1, silo=0):
	"""
		Parser for oct2grid
	"""

	N = pow(2,level)
	if xmax == -1 :
		xmax = N
	if ymax == -1 :
		ymax = N
	if zmax == -1 :
		zmax = N
		

	f = field2ID(field)

	print "********** Oct2grid **********"
	print "File = ", data
	print "Field = ", field
	print "Level = ", level
	print "******************************"

	
	dx =pow(2.,-level)  
	
	xmin *= dx
	xmax *= dx	
	
	ymin *= dx
	ymax *= dx	
	
	zmin *= dx
	zmax *= dx
		
	o2g  = "utils/oct2grid "
	
	out  = data.replace("grid","cube"+ str(level)) + "." + field
	commande = (o2g + " " + 
							data + " " +  
							str(level) + " "  + 
							str(f) + " " + 
							out + " " + 
							str(nproc) + " " + 
							str(silo) + 
							" -1 " + 
							str(xmin) + " " + 
							str(xmax) + " " + 
							str(ymin) + " " + 
							str(ymax) + " " + 
							str(zmin) + " " + 
							str(zmax)  )

	if os.path.isfile(out):
		print "oct2grid : file allready exist"
		if force:
			print "but do it anyway!"
			os.system(commande)
	else :
		os.system(commande)

