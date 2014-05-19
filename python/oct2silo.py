#!/usr/bin/env python
from fonction_IO  import *
from fonction_amr import *


if __name__ == "__main__":	


	args = getargs()
	args.field = ["rfield.temp"]
	filename = args.files[0][:-7]
	filename = filename[:-10] +  "grid" +  filename[-6:] 
	denoct2grid(filename, args, 1)
