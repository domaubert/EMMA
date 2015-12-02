#------------ MAIN OPTIONS --------------------
DEFINES  +=  -DPIC
DEFINES  +=  -DWHYDRO2
DEFINES  +=  -DWGRAV
DEFINES  +=  -DWRAD
DEFINES  +=  -DWRADHYD
DEFINES  +=  -DTESTCOSMO
DEFINES  +=  -DSTARS
DEFINES  +=  -DSUPERNOVAE
#DEFINES  +=  -DZOOM

#------------ PRECISION OPTIONS ---------------------
#DEFINES  +=  -DSINGLEPRECISION

#------------ MPI OPTIONS ---------------------
DEFINES  +=  -DWMPI
#DEFINES  +=  -DFLOORDT

#------------ OPEN MP OPTIONS ---------------------
#DEFINES += -DWOMP

#------------ CUDA OPTIONS ---------------------
DEFINES  +=  -DWCUDA_ERR
#DEFINES  +=  -DNOCOMP

#------------ ICs OPTIONS ---------------------
#DEFINES  +=  -DJUSTIC
#DEFINES  +=  -DGRAFIC -DBULKFLOW
DEFINES  +=  -DGRAFIC
#DEFINES  +=  -DZELDOVICH
#DEFINES  +=  -DEVRARD
#DEFINES  +=  -DEDBERT
#DEFINES  +=  -DTUBE
#DEFINES  +=  -DPARTN
#DEFINES  +=  -DPART2
#DEFINES  +=  -DWRADTEST
#DEFINES  +=  -DTESTCLUMP # RADTEST MUST BE SET
#DEFINES  +=  -DSNTEST    # RADTEST MUST BE SET

#------------ PIC OPTIONS ----------------------
#DEFINES += -DPART_EGY
#DEFINES += -DPERFECT
#DEFINES += -DAMRPART

#------------ GRAV OPTIONS ----------------------
#DEFINES  +=  -DFASTGRAV
DEFINES  += -DONFLYRED

# ----------- HYDRODYNAMICS OPTIONS ------------
DEFINES  +=  -DRIEMANN_HLLC
#DEFINES  +=  -DRIEMANN_EXACT
DEFINES  +=  -DPRIMITIVE
DEFINES  +=  -DDUAL_E
#DEFINES  +=  -DNOADX

# ----------- RADIATION OPTIONS ------------
DEFINES  += -DWCHEM
#DEFINES  += -DHELIUM
DEFINES  += -DS_100000
#DEFINES  += -DSTARBURST_5MYR
DEFINES  += -DCOOLING
#DEFINES  += -DUVBKG
DEFINES  += -DSEMI_IMPLICIT
DEFINES  += -DACCEL_RAD_STAR
DEFINES += -DOTSA
#DEFINES  += -DHOMOSOURCE
#DEFINES  += -DRADSTEP
DEFINES  += -DCOARSERAD

# ----------- STAR FORMATION OPTIONS ------------
#DEFINES  += -DSIMPLESTAR // simple hack to perform SF � la Rasera
#DEFINES  += -DJEANSCRIT // adds Jeans criterion to trigger SF
#DEFINES  += -DSCHAYE // use Schaye Polytropic law and SF laws
DEFINES  += -DDECREASE_EMMISIVITY_AFTER_TLIFE
#DEFINES  += -DMULTIPLESTARS // N part of mass M instead of 1 part of mass N*M
#DEFINES  += -DCONTINUOUSSTARS // All cells allowed to form star, effectivelly form stars  
#DEFINES  += -DTEMPCRIT // add a 2e3K temperature criterion

# ---- BOUNDARY CONDITIONS (PERIODIC BY DEFAULT)--
#DEFINES  +=  -DTRANSZM
#DEFINES  +=  -DTRANSZP
#DEFINES  +=  -DTRANSYM
#DEFINES  +=  -DTRANSYP
#DEFINES  +=  -DTRANSXM
#DEFINES  +=  -DTRANSXP
#DEFINES  +=  -DREFXM # TRANS must be turned on too
#DEFINES  +=  -DREFYM # TRANS must be turned on too
#DEFINES  +=  -DREFZM # TRANS must be turned on too

# ---- OUTPUT--------------
#DEFINES += -DWDBG
#DEFINES += -DMOVIE
DEFINES += -DBKP
DEFINES += -DMULTIFOLDER
#DEFINES += -DALLOCT
#DEFINES += -DMPIIO
