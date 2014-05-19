#!/usr/bin/env python

if __name__ == "__main__":


	folder= "./"
	filename = folder + "star"



	try:
		f = open(filename, "r")	
	except IOError:
		print 'cannot open', filename 

	data = f.read()
	f.close()
	print "read OK"
	

	data = [x for x in data.split('\n')]
	data = data[2:]

	for i in range(len(data)):
		data[i] = data[i][8:]
		data[i] = [x for x in data[i].split('  ')] 
		if (len(data[i]) != 7 ) :
			data[i] = 0


############################################################################
	
	filename = filename + ".3D"
	try:
		f = open(filename, "wb")	
		print "writing", filename

	except IOError:
		print 'cannot open', filename 

	f.write("x y z value\n");

	inc = 1
	p=0
	while p<len(data):
		if data[p] : 
			f.write("%s %s %s %s\n" % (data[p][2],data[p][3],data[p][4],data[p][5]))
		p += inc

	f.close()	





