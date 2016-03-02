!------------------------------------------------------------------------
! This parallel Friend of Friend implementation is based on a sequential
! implementation written by Edouard Audit (CEA).
! This parallel implementation has been written by Fabrice Roy
! CNRS / LUTH (Observatoire de Paris)
! mail: fabrice.roy@obspm.fr
!------------------------------------------------------------------------


Program friend

  Use modconst
  Use modvariable
  Use modparam
  Use modmpicom
  Use modfofpara
  Use modio
  Use modtiming

  Implicit none

  ! Local variables
  
  Integer :: errcode
  Character(len=90):: filelog
  
  !----------------------------------------------------

  ! Initialization of MPI
  Call Mpi_Init(mpierr)
  Call Mpi_Comm_size(Mpi_Comm_World, procNB, mpierr)

  ! Creation of a grid of processes.
  ! The grid is a cube, there are procNB**(1/3) processes in each dimension.
  ! The grid is periodic in each dimension.
  dims = int(procNB**(1./3.))
  periods = .true. 
  
  Call Mpi_Cart_create(Mpi_Comm_World,3,dims,periods,.true.,MPICube,mpierr)

  ! Process ID and number 
  Call Mpi_Comm_rank  (MPICube, procID, mpierr)
  Call Mpi_Cart_coords(MPICube, procID, 3, CubeCoord, mpierr)

  ! Definition des voisins pour les 6 directions avec des conditions au bord periodiques
  ! Intialization of the 'neighboor processes id' 
  ! Neighboor : 1 = backward
  !             2 = forward
  !             3 = left
  !             4 = right
  !             5 = down
  !             6 = up

  Call Mpi_Cart_shift(MPICube,0,1,voisin(1),voisin(2),mpierr)
  Call Mpi_Cart_shift(MPICube,1,1,voisin(3),voisin(4),mpierr)
  Call Mpi_Cart_shift(MPICube,2,1,voisin(5),voisin(6),mpierr)


  ! Memory allocation for the input parameters broadcast buffer
  ! 90*4               : 4x Character(len=90) root, pathinput, namepart, nameinfo 
  ! 3                  : 1x Character(len=3)  code_index
  ! 4*bit_size(Mmin)/8 : 3x Integer(kind=4)   grpsize, snapnum,Mmin, Mmax
  ! kind(perco)        : 1x Real(kind=4)      perco
  ! 6*4                : 6x Logical           star, metal, outcube, dofof, readfromcube, dotimings

  h_length = 90*4 + 3 + 4* bit_size(Mmin)/8 + kind(perco) + 6*4 !HACKDOM
  allocate (header(0:h_length-1))
  h_pos = 0


  ! The process 0 read the input parameters and pack them for the broadcast.
  If(procID==0) Then
     
     ! Read input parameters in file 'fof.in'
     Open(10,file='fof.in',iostat=errcode)
     If(errcode>0) Then
        Print *,'** Error opening input file fof.in. Please check this file. **'
        Call Mpi_Abort(Mpi_Comm_World,errcode,mpierr)
     End If
     
     Read(10,'(A90)') root
     Read(10,*) code_index
     Read(10,'(A90)') pathinput
     Read(10,'(A90)') namepart
     Read(10,'(A90)') nameinfo
     Read(10,*) grpsize
     Read(10,*) perco
     Read(10,*) snapnum ! HACK DOM
     Read(10,*) Mmin
     Read(10,*) Mmax
     Read(10,*) star
     Read(10,*) metal
     Read(10,*) outcube
     Read(10,*) dofof
     Read(10,*) readfromcube
     Read(10,*) dotimings
     Close(10)


     Print *,'Parallel FoF'
     Print *,procNB,' processes:'

     ! Open log file
     filelog = trim(root)//'.log'
     Open(Unit=Ulog,file=filelog)
     Write(Ulog,*) 'Parallel FoF'
     Write(Ulog,*) procNB,' processes:'

     ! Write input parameters in a .inp file
     Call outputparameters()

     ! Pack input parameters in the buffer
     Call Mpi_Pack(        root, 90, Mpi_Character, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(  code_index,  3, Mpi_Character, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(   pathinput, 90, Mpi_Character, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(    namepart, 90, Mpi_Character, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(    nameinfo, 90, Mpi_Character, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(     grpsize,  1,   Mpi_Integer, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(       perco,  1,      Mpi_Real, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(     snapnum,  1,   Mpi_Integer, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(        Mmin,  1,   Mpi_Integer, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(        Mmax,  1,   Mpi_Integer, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(        star,  1,   Mpi_Logical, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(       metal,  1,   Mpi_Logical, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(     outcube,  1,   Mpi_Logical, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(       dofof,  1,   Mpi_Logical, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(readfromcube,  1,   Mpi_Logical, header, h_length, h_pos, Mpi_Comm_World, mpierr)
     Call Mpi_Pack(   dotimings,  1,   Mpi_Logical, header, h_length, h_pos, Mpi_Comm_World, mpierr)

     ! Print input parameters on screen
     Print *, ''
     Print *, "Root:                                           ",trim(root)
     Print *, "Type of input files:                            ",code_index
     Print *, "Path to input files:                            ",trim(pathinput)
     Print *, "Particle files base name:                       ",trim(namepart)
     Print *, "Info file base name:                            ",trim(nameinfo)
     Print *, "Size of groups of inputs:                       ",grpsize
     Print *, "Percolation parameter:                          ",perco
     Print *, "Snapshot number      :                          ",snapnum
     Print *, "Minimum mass of halo to be analyzed:            ",Mmin
     Print *, "Maximum mass of halo to be analyzed:            ",Mmax
     Print *, "Were stars written in RAMSES files:             ",star
     Print *, "Were metallicities written in RAMSES files:     ",metal
     Print *, "Write cubes of particles:                       ",outcube
     Print *, "Perform friends of friends halo detection:      ",dofof
     Print *, "Read particles from cube files:                 ",readfromcube
     Print *, "Perform timings (imply extra synchronisations): ",dotimings
     Print *, ''

     ! Write input parameters in .log file
     Write(Ulog,*) 'INPUT PARAMETERS fof.in'
     Write(Ulog,*) root
     Write(Ulog,*) code_index
     Write(Ulog,*) pathinput
     Write(Ulog,*) namepart
     Write(Ulog,*) nameinfo
     Write(Ulog,*) grpsize
     Write(Ulog,*) perco
     Write(Ulog,*) snapnum
     Write(Ulog,*) Mmin
     Write(Ulog,*) Mmax
     Write(Ulog,*) star
     Write(Ulog,*) metal
     Write(Ulog,*) outcube
     Write(Ulog,*) dofof
     Write(Ulog,*) readfromcube
     Write(Ulog,*) dotimings

  End If

  ! Broadcast of input parameters
  Call Mpi_Bcast(header,h_length,Mpi_Packed,0,Mpi_Comm_World,mpierr)

  ! Processes with ID != 0 unpack the input parameters
  If(procID /= 0) Then
     Call Mpi_Unpack(header, h_length, h_pos,         root, 90, Mpi_Character, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,   code_index,  3, Mpi_Character, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,    pathinput, 90, Mpi_Character, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,     namepart, 90, Mpi_Character, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,     nameinfo, 90, Mpi_Character, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,      grpsize,  1,   Mpi_Integer, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,        perco,  1,      Mpi_Real, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,      snapnum,  1,   Mpi_Integer, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,         Mmin,  1,   Mpi_Integer, Mpi_Comm_World, mpierr)
     Call Mpi_UnPack(header, h_length, h_pos,         Mmax,  1,   Mpi_Integer, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,         star,  1,   Mpi_Logical, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,        metal,  1,   Mpi_Logical, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,      outcube,  1,   Mpi_Logical, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,        dofof,  1,   Mpi_Logical, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos, readfromcube,  1,   Mpi_Logical, Mpi_Comm_World, mpierr)
     Call Mpi_Unpack(header, h_length, h_pos,    dotimings,  1,   Mpi_Logical, Mpi_Comm_World, mpierr)
  End If

  ! broadcast buffer deallocated
  Deallocate(header)

  ! code_index determines the type of cosmological simulation to be analyzed:
  ! RA2 = RAMSES v2
  ! RA3 = RAMSES v3
  ! Read the output of the cosmological simulation
  If (code_index.eq.'RA2' .or. code_index.eq.'RA3'.or. code_index.eq.'EMM') Then
     !Call RAMSES_lecture()
     call EMMA_lecture()
  Else
     If(procID==0) Print *,'** Wrong file type. Possibilities are RA2 or RA3. OR EMMA **'
     errcode = 2
     Call Mpi_Abort(Mpi_Comm_World,errcode,mpierr)
     Stop
  End If

  ! Print on screen and in log file 
  If(procID==0) Then
     Write(* ,*)   'nz = ',nres,'ngrid = ',ngrid
     Write(Ulog,*) 'nz = ',nres,'ngrid = ',ngrid
  End If

  ! Write cube of particles files if requested
  If(outcube) Then
     If(procID==0) Print *,      'Write particles distributed in a cartesian grid'
     If(procID==0) Write(Ulog,*) 'Write particles distributed in a cartesian grid'
     Call outputcube()
  End If

  
  If(procID==0) Then
     Print *,' '
     Print *,'Friends of Friends halo detection'
     Print *,' '
     Write(Ulog,*) 'Friends of Friends halo detection'
  End If


  ! Parallele Friends of Friends if requested
  If(dofof) Then
     Call parafof()
  Else
     tFoF     = 0.0
     tFoFinit = 0.0
     tFoFloc  = 0.0
     tRaccord = 0.0
     tObs     = 0.0
     tOut     = 0.0
  End If

  ! Output of timings if requested
  If(dotimings .and. procID==0) Then
     Print *,'Friend of Friend termine'
     Print *,''
     Print *,'Temps d''execution:'
     Print *,'Lecture:',tReadRA,' dont'
     Print *,'        initialisation        :',tInitRead
     Print *,'        lecture des fichiers  :',tReadFile
     Print *,'        partage des particules:',tTailPart
     Print *,''
     Print *,'Friend of Friend:',tFoF,'dont'
     Print *,'        initialisatin:',tFoFinit
     Print *,'        FoF local    :',tFoFloc
     Print *,'        raccordement :',tRaccord
     Print *,'        calcul d''observables:',tObs
     Print *,'        sorties:',tOut

     Write(Ulog,*) 'Friend of Friend termine'
     Write(Ulog,*) ''
     Write(Ulog,*) 'Temps d''execution:'
     Write(Ulog,*) 'Lecture:',tReadRA,' dont'
     Write(Ulog,*) '        initialisation        :',tInitRead
     Write(Ulog,*) '        lecture des fichiers  :',tReadFile
     Write(Ulog,*) '        partage des particules:',tTailPart
     Write(Ulog,*) ''
     Write(Ulog,*) 'Friend of Friend:',tFoF,'dont'
     Write(Ulog,*) '        initialisatin:',tFoFinit
     Write(Ulog,*) '        FoF local    :',tFoFloc
     Write(Ulog,*) '        raccordement :',tRaccord
     Write(Ulog,*) '        calcul d''observables:',tObs
     Write(Ulog,*) '        sorties:',tOut
  End If

  ! Deallocation of x,v,id which are allocated in RAMSES_lecture subroutine
  Deallocate(x,v,id)

  ! Close log file
  Close(Ulog)

  ! Finalize MPI
  Call Mpi_Finalize(mpierr)

End Program friend
