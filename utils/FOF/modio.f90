Module modio

    Use modconst

    Integer(kind=4)   :: mynpart       ! local particle number
    Integer(kind=4)   :: ndim          ! number of dimensions
    Integer(kind=4)   :: lmin          ! minimum mesh refinement level
    Integer(kind=4)   :: nres          ! 1-D resolution: number of grid points in each dimension ; nres = 2 ** lmin
    Integer(kind=PRI) :: nptot         ! total particle number: nptot = nres ** 3
    Integer(kind=PRI) :: ngrid         ! total number of grid points: ngrid = nres ** 3
    Real   (kind=SP)  :: xmin, xmax, ymin, ymax, zmin, zmax  ! min and max (x,y,z) for each process
    Integer(kind=PRI) :: nptot_grp_threshold=10000000000 !!!!RY group modif nptot threshold for fof for writing in grp.
  Contains


    ! Write input parameters in file .inp
    Subroutine outputparameters()
        Use modparam
        Use modmpicom
        Implicit none

        Character(len=90) :: fileopa

        ! Open file .inp and write input parameters
        fileopa = trim(root)//'.inp'
        Open(Unit=Uopa,file=fileopa)
        Write(Uopa,*) 'Nb of processes:'
        Write(Uopa,*) procNB
        Write(Uopa,*) ''
        Write(Uopa,*) 'Input parameters:'
        Write(Uopa,*) root
        Write(Uopa,*) code_index
        Write(Uopa,*) pathinput
        Write(Uopa,*) namepart
        Write(Uopa,*) nameinfo
        Write(Uopa,*) grpsize
        Write(Uopa,*) perco
        Write(Uopa,*) Mmin
        Write(Uopa,*) Mmax
        Write(Uopa,*) star
        Write(Uopa,*) metal
        Write(Uopa,*) outcube
        Write(Uopa,*) dofof
        Write(Uopa,*) readfromcube
        Write(Uopa,*) dotimings

        Close(Uopa)

    End Subroutine outputparameters


    ! Read particles files created by RAMSES
    ! RAMSES writes as many particles files as it used process to run the simulation.
    ! Each process will read a part of the files.
    ! Particles are spread among the different files and have to be sorted by position and distributed between the processes.
    ! Each process in the parallel FoF will analize a cube of particles, this means the particles have to be geographically
    ! distributed.
    Subroutine ramses_lecture()
        Use modvariable
        Use modparam
        Use modmpicom
        Use modtiming
        Implicit none

        ! Local variables
        Character(len=5)               :: ncharcpu
        Character(len=9)               :: tmpstr1, tmpstr2
        Character(len=200)             :: nomfich
        Character(len=13)              :: dumchar
        Character(len=11)              :: grpchar

        Integer(kind=4)                :: i, j, icpu, idim   ! loop variables
        Integer(kind=4)                :: destCoord(3)       ! coords of the destination MPI process in MPI process cart
        Integer(kind=4)                :: nrecv              ! number of elements received in a Mpi_Recv
        Integer(kind=4)                :: recvpoint          ! address of the 1st element received in the local vector
        Integer(kind=4)                :: mynbfile            ! number of RAMSES part files read by local process
        Integer(kind=4)                :: firstp, lastp      ! id of 1st and last RAMSES part file read
        Integer(kind=4), allocatable   :: npartvloc(:), npartv(:)  ! temp and global table of particle numbers for each process
        Integer(kind=4)                :: n_i, n_j, n_k, nsd, ind
        Integer(kind=4)                :: nproc       ! process number  read in RAMSES info file
        Integer(kind=4)                :: ncpu2       ! process number  read in RAMSES part files
        Integer(kind=4)                :: ndim2       ! dimension       read in RAMSES part files
        Integer(kind=4)                :: npartloc    ! particle number read in RAMSES part files
        Integer(kind=4)                :: prov, dest  ! provenance and destination process number for p2p MPI communications
        Integer(kind=4)                :: mpistat(Mpi_Status_Size)   ! status of MPI communication
        Integer(kind=PRI)              :: tmplongint              ! temp integer8 variable
        Integer(kind=PRI), allocatable :: tmpi(:), tmpsendi(:)    ! TYPE VARIABLE EN FONCTION DU NB DE PART
        Integer(kind=4)                :: errcode
        Integer(kind=4)                :: grpnb

        Real(kind=SP), allocatable     :: tmpsimple(:)            ! temporary variable for Ramses v2 output
        Real(kind=DP), allocatable     :: tmpdouble(:)            ! temporary variable for Ramses v3 output
        Real(kind=SP), allocatable     :: tmpsendx(:,:),tmpsendv(:,:)
        Real(kind=SP), allocatable     :: tmpx(:,:), tmpv(:,:)
        Real(kind=SP)                  :: deltasd

        integer(kind=4)                ::nmod
        
        ! buffer for Mpi_Pack
        ! 3*bit_size(nres)/8 : 3x Integer(kind=4) nproc, ndim, nres
        h_length = 3* bit_size(nres)/8
        Allocate (header(0:h_length-1))
        h_pos = 0

        ! Timer initialization
        time0 = Mpi_Wtime()

        ! group directory base name: particles files are dispatched in different 'group' directories
        grpchar = 'group_00001'

        ! Process 0 reads parameters in info file
        If(procID == 0) Then
            ! RAMSES v2 particles files
            If(code_index.eq.'RA2') Then
                Print *,'Reading Ramses v2 output...'
                Write(Ulog,*) 'Reading Ramses v2 output...'
                ! RAMSES v3 particles files
            Else if(code_index.eq.'RA3') Then
                Print *,'Reading Ramses v3 output...'
                Write(Ulog,*) 'Reading Ramses v3 output...'
            End If

            ! info file
            !!!!RY group modif
            !nomfich = trim(pathinput)//'/'//trim(grpchar)//'/'//trim(nameinfo)
            if(grpsize>0) then
               nomfich = trim(pathinput)//'/'//trim(grpchar)//'/'//trim(nameinfo)
            else
               nomfich = trim(pathinput)//'/'//trim(nameinfo)
            endif
            !!!END RY group modif
            Print *,'Reading RAMSES info file:',trim(nomfich)

            Open(11,file=nomfich, form='formatted', Status='Old', Iostat=errcode)
            If(errcode > 0) Then
                Call EmergencyStop('Error opening '//trim(nameinfo)//' file',errcode)
            End If
            Rewind 11

            ! Read number of particles files, nb of dimensions and level of refinement in the coarse grid (RAMSES is based on AMR)
            Read(11,'(A13,I11)') dumchar,nproc
            Read(11,'(A13,I11)') dumchar,ndim
            Read(11,'(A13,I11)') dumchar,lmin

            Close(11)

            ! nb of grid point in each dimension = 2^lmin
            nres = 2**lmin

            Write(*,*) 'Number of:'
            Write(*,'(A25,I6)') ' - files for each output:',nproc
            Write(*,'(A25,I6)') ' - dimensions:           ',ndim
            Write(*,'(A25,I6)') ' - grid points:          ',nres
            Write(Ulog,*) 'nb_proc = ',nproc,'ndim = ',ndim,'nres = ',nres

            ! Pack these 3 parameters for broadcast
            Call Mpi_Pack(nproc, 1, Mpi_Integer, header, h_length, h_pos, Mpi_Comm_World, mpierr)
            Call Mpi_Pack( ndim, 1, Mpi_Integer, header, h_length, h_pos, Mpi_Comm_World, mpierr)
            Call Mpi_Pack( nres, 1, Mpi_Integer, header, h_length, h_pos, Mpi_Comm_World, mpierr)

        End If

        ! Broadcast
        Call Mpi_Bcast(header,h_length,Mpi_Packed,0,Mpi_Comm_World,mpierr)

        ! Other processes unpack the 3 parameters
        If(procID /= 0) Then

            Call Mpi_Unpack(header, h_length, h_pos, nproc, 1, Mpi_Integer, Mpi_Comm_World, mpierr)
            Call Mpi_Unpack(header, h_length, h_pos,  ndim, 1, Mpi_Integer, Mpi_Comm_World, mpierr)
            Call Mpi_Unpack(header, h_length, h_pos,  nres, 1, Mpi_Integer, Mpi_Comm_World, mpierr)

        End If

        ! total nb of grid points
        ngrid = int(nres,kind=8)**3

        If(allocated(header)) Deallocate(header)

        If(procID==0) Print *,'Reading positions...'
        If(procID==0) Write(Ulog,*) 'Reading positions...'

        ! total nb of particles is initialized to 0
        nptot = 0
        ! each process read mynbfile particles files
        mynbfile = nproc / procNB
        ! index of first and last particles file read by current process
        firstp  = procID * mynbfile + 1
        lastp   = (procID+1) * mynbfile
        
        !!!!When mod(nproc, procNB) .ne. 0!!!RY
        nmod=mod(nproc,procNB)
        if (nmod.ne.0)then
           if(procID.le.nmod-1)then
              mynbfile=mynbfile+1
              firstp  = procID * mynbfile + 1
              lastp   = (procID+1) * mynbfile
           else
              firstp  = procID * mynbfile + 1+nmod
              lastp   = (procID+1) * mynbfile+nmod
           endif
        endif
        !!!!End RY!!!

        ! nb of particles in files read by current process
        mynpart = 0

        ! particles will  be distributed "geographicaly" among the process, i.e. at the end of this subroutine
        ! each process will keep only the particles located in a particular subdomain of the whole simulation domain
        ! allocate array containing the nb of particles that will remain on each process at the end of this subroutine (npartv)
        ! and array containing the nb of particles destined to each process but read by current process (npartvloc)
        Allocate(npartv(procNB))
        Allocate(npartvloc(procNB))

        npartv = 0
        npartvloc = 0

        ! nsd = nb of subdomains in each direction, i.e. there are nsd**3 subdomains
        nsd = int(procNB**(1./3.))
        ! deltasd = dimension of a subdomain, i.e. each subdomain is a cube of size deltasd
        deltasd = 1./nsd

        If(procID == 0) Then
            Write(*,*) 'Number of subdomains in each dimension:',nsd
            Write(*,*) 'Size of each subdomain:',deltasd
        End If

        ! definition of the current subdomain by min and max coordinates
        xmin =  CubeCoord(1)      * deltasd
        xmax = (CubeCoord(1) + 1) * deltasd
        ymin =  CubeCoord(2)      * deltasd
        ymax = (CubeCoord(2) + 1) * deltasd
        zmin =  CubeCoord(3)      * deltasd
        zmax = (CubeCoord(3) + 1) * deltasd

        ! loop over the "file/cpu" index characterizing the input particles files to be read by current process
        Do icpu = firstp,lastp
            Write(ncharcpu(1:5),'(I5.5)') icpu
            ! the RAMSES files are written in group_XXXXX directories, where XXXXX is the group number
            ! the grpsize first files are in group number 1, the grpsize following are in group number 2, etc...
            
            !!!!RY group modif
            !grpnb = (icpu-1)/grpsize + 1
            !Write(grpchar(7:11),'(I5.5)') grpnb
            !nomfich = trim(pathinput)//'/'//trim(grpchar)//'/'//trim(namepart)//trim(ncharcpu)
            if(grpsize >0)then
               grpnb = (icpu-1)/grpsize + 1
               Write(grpchar(7:11),'(I5.5)') grpnb
               nomfich = trim(pathinput)//'/'//trim(grpchar)//'/'//trim(namepart)//trim(ncharcpu)
            else
               nomfich = trim(pathinput)//'/'//trim(namepart)//trim(ncharcpu)
            endif
            !!!!END RY group modif
            ! read number of files, dimensions and particles
            Open(unit=1,file=nomfich,status='old',form='unformatted',Iostat=errcode)
            If(errcode > 0) Then
               Call EmergencyStop('Error opening '//trim(nameinfo)//' file',errcode)
            End If
            Read(1,Iostat=errcode) ncpu2
            If(errcode/=0) Print *,'Erreur de lecture: ',errcode
            Read(1,Iostat=errcode) ndim2
            If(errcode/=0) Print *,'Erreur de lecture: ',errcode
            Read(1,Iostat=errcode) npartloc
            If(errcode/=0) Print *,'Erreur de lecture: ',errcode
            Close(1)

            ! check if nb of files and dimensions are consistent with what is written in info file
            If((ncpu2/=nproc).Or.(ndim2/=ndim)) Then
                Call EmergencyStop('Files'//trim(nomfich)// ' and '//trim(nameinfo)//' are not consistent for ncpu and/or ndim',2)
            End If

            ! nb of particles read by current process
            mynpart = mynpart + npartloc
        End Do

        ! we use a temporary long integer in case of simulations larger than 1024**3
        tmplongint = mynpart
        ! each mynpart is added -the resulst is nptot, total nb of particles - and broadcasted in a single collective communication
        Call Mpi_AllReduce(tmplongint,nptot,1,MPI_PRI,Mpi_Sum,Mpi_Comm_World,mpierr)

        If(procID == 0) Then
            Write(* ,*)'There are ',nptot,' DM particles'
            Write(Ulog,*)'There are ',nptot,' DM particles'
        End If

        ! memory allocation for temporary arrays which will contain positions, velocities and id's
        ! of the particles read by current process
        Allocate(tmpx(3,mynpart))
        Allocate(tmpv(3,mynpart))
        Allocate(tmpi(mynpart))

        ! intialization
        tmpx=0.
        mynpart = 0

        ! Begin timings if requested
        If(dotimings) Then
            Call Mpi_Barrier(MPICube,mpierr)
            timeInt = Mpi_Wtime()
            tInitRead = timeInt - time0
        End If

        ! current process reads its particles files again, but this time it reads position, velocities and id's for each particle
        ! in the file
        Do icpu = firstp,lastp
           !!!!RY group modif
           ! Write(ncharcpu(1:5),'(I5.5)') icpu
           ! grpnb = (icpu-1)/grpsize + 1
           ! Write(grpchar(7:11),'(I5.5)') grpnb
           ! nomfich = trim(pathinput)//'/'//trim(grpchar)//'/'//trim(namepart)//trim(ncharcpu)
           Write(ncharcpu(1:5),'(I5.5)') icpu
            if(grpsize >0) then
               grpnb = (icpu-1)/grpsize + 1
               Write(grpchar(7:11),'(I5.5)') grpnb
               nomfich = trim(pathinput)//'/'//trim(grpchar)//'/'//trim(namepart)//trim(ncharcpu)
            else
               nomfich = trim(pathinput)//'/'//trim(namepart)//trim(ncharcpu)
            endif
            !!!END RY group modif
           
            Open(unit=1,file=nomfich,status='old',form='unformatted')
            Read(1) ncpu2
            Read(1) ndim2
            Read(1) npartloc

            ! There are some differences between RAMSES v2 and RAMSES v3 files
            ramsesV2 : If(code_index.eq.'RA2') Then
            ! allocate an array of real simple precision to read position/velocities
            Allocate(tmpsimple(1:npartloc))
            ! Read positions
            Do idim = 1,ndim
                Read(1) tmpsimple
                ! put all positions in tmpx vector
                tmpx(idim,mynpart+1:mynpart+npartloc) = tmpsimple
            End Do

            ! Read velocities in a dummy variable
            Do idim = 1,ndim
                Read(1) tmpsimple
                ! put all velocities in tmpv vector
                tmpv(idim,mynpart+1:mynpart+npartloc) = tmpsimple
            End Do

            ! Read masses in a dummy variable
            Read(1) tmpsimple
            Deallocate(tmpsimple)

            End If ramsesV2

            ! There are some differences between RAMSES v2 and RAMSES v3 files
            ramsesV3 : If(code_index.eq.'RA3') Then
            ! allocate an array of real double precision to read position/velocities
            ! these double precision values will be stored in simple precision variables
            Allocate(tmpdouble(1:npartloc))
            ! Read some dummy variables
            Read(1)
            Read(1)
            Read(1)
            Read(1)
            Read(1)

            ! Read positions
            Do idim = 1,ndim
                Read(1) tmpdouble
                ! put all positions in tmpx vector
                tmpx(idim,mynpart+1:mynpart+npartloc) = tmpdouble
            End Do

            ! Read velocities in a dummy variable
            Do idim = 1,ndim
                Read(1) tmpdouble
                ! put all velocities in tmpv vector
                tmpv(idim,mynpart+1:mynpart+npartloc) = tmpdouble
            End Do

            ! Read masses in a dummy variable
            Read(1) tmpdouble

            If(star) Then
                Read(1) tmpdouble
                If(metal) Read(1) tmpdouble
            End If

            Deallocate(tmpdouble)
            End If ramsesV3
            ! End of the differences between the two versions of RAMSES

            ! Read particle id
            Read(1) tmpi(mynpart+1:mynpart+npartloc)            

            ! loop over the particles read in the current particles file
            Do j = mynpart+1,mynpart+npartloc
                ! particle positions must be >= 0 and <1.0
                ! as the domain is periodic, every position=1.0 is set to 0.0
                Do idim = 1,ndim
                    If(tmpx(idim,j)==1.0) tmpx(idim,j) = 0.0
                End Do
                ! in which subdomain is located this particle?
                n_i = int(tmpx(1,j)/deltasd)
                n_j = int(tmpx(2,j)/deltasd)
                n_k = int(tmpx(3,j)/deltasd)
                ind = nsd**2 *n_i + nsd*n_j + n_k + 1
                ! this means subdomain number 'ind' has one more particle to take care of
                npartvloc(ind) = npartvloc(ind)+1
            End Do

            ! close particles file
            Close(1)
            ! current process has read mynpart particles
            mynpart = mynpart+npartloc
        End Do

        ! add and broadcast npartvloc so that each process knows how many particles it and its friends will have to analize
        Call Mpi_AllReduce(npartvloc,npartv,procNB,Mpi_Integer,Mpi_Sum,Mpi_Comm_World,mpierr)

        ! timing if requested
        If(dotimings) Then
            tReadfile = Mpi_Wtime() - timeInt
            timeInt = Mpi_Wtime()
        End If

        ! ------------------------------------------------
        ! sharing out of the particles between processes
        ! ------------------------------------------------

        ! allocate final positions/velocities/id's arrays
        Allocate (x(3,npartv(procID+1)))
        Allocate (v(3,npartv(procID+1)))
        Allocate (id(npartv(procID+1)))

        ! index of the position in the arrays where the data received from other processes must be stored
        recvpoint = 1

        ! loop on processes: each process will received from the other processes the positions/velocities/id's of the particles
        ! it has to analyze
        processus : Do i = 1,procNB - 1
        ! ID of the process current process has to send data to
        dest = mod(procID + i,procNB)
        ! ID of the process current process has to receive data from
        prov = mod(procID + procNB - i, procNB)
        ! Coordinates of dest process in the process grid
        Call Mpi_Cart_coords(MPICube,dest,3,destCoord,mpierr)
        ! current process sends and receives the size of the data it will send/receive
        Call Mpi_Isend(npartvloc(dest+1),1,Mpi_Integer,dest,procID,MpiCube,mpireqs1,mpierr)
        Call Mpi_Irecv(nrecv,1,Mpi_Integer,prov,prov,MpiCube,mpireqr1,mpierr)
        ! current process determines the boundaries of the region it has to send to dest
        xmin =  destCoord(1)      * deltasd
        xmax = (destCoord(1) + 1) * deltasd
        ymin =  destCoord(2)      * deltasd
        ymax = (destCoord(2) + 1) * deltasd
        zmin =  destCoord(3)      * deltasd
        zmax = (destCoord(3) + 1) * deltasd

        ! there are npartvloc(dest+1) particles on the current process that it has to send to dest
        ! if this nb is  greater than 0 then current process proceeds
        If(npartvloc(dest+1)>0) Then
            ! allocation of temporary arrays
            Allocate(tmpsendx(3,npartvloc(dest+1)))
            Allocate(tmpsendv(3,npartvloc(dest+1)))
            Allocate(tmpsendi(  npartvloc(dest+1)))

            ind = 1
            ! loop over local particles
            Do j=1,mynpart
                ! if the particle is located in the cube of dest, then it is added to the send buffers
                If(tmpx(1,j)>=xmin .and. tmpx(1,j)<xmax .and. &
                tmpx(2,j)>= ymin .and. tmpx(2,j) < ymax .and. &
                tmpx(3,j)>= zmin .and. tmpx(3,j) < zmax) Then
                tmpsendx(:,ind) = tmpx(:,j)
                tmpsendv(:,ind) = tmpv(:,j)
                tmpsendi(  ind) = tmpi(j)
                ind=ind+1
                End If
            End Do

            ! if current process has not found the same nb of particles in dest cube as it had found while
            ! reading from files, then there is a bug somewhere
            If(ind/=npartvloc(dest+1)+1) Then
                Print *,'Erreur dans la repartition des particules'
                Call Mpi_Finalize(mpierr)
                Stop
            End If

       End If

       ! wait for size of send/recv arrays communication to complete
            Call Mpi_Wait(mpireqs1,mpistat,mpierr)
            Call Mpi_Wait(mpireqr1,mpistat,mpierr)

            ! send particles to dest process if needed
            If(npartvloc(dest+1)/=0) Then
                Call Mpi_Isend(tmpsendx,3*npartvloc(dest+1),Mpi_Real,dest,procID,MpiCube,mpireqs1,mpierr)
                Call Mpi_Isend(tmpsendv,3*npartvloc(dest+1),Mpi_Real,dest,  dest,MpiCube,mpireqs2,mpierr)
                Call Mpi_Isend(tmpsendi,  npartvloc(dest+1), MPI_PRI,dest,procID,MpiCube,mpireqs3,mpierr)
            End If
            ! receive particles from recv process if needed
            If(nrecv/=0) Then
                Call Mpi_Irecv(x(1,recvpoint),3*nrecv,Mpi_Real,prov,  prov,MpiCube,mpireqr1,mpierr)
                Call Mpi_Irecv(v(1,recvpoint),3*nrecv,Mpi_Real,prov,procID,MpiCube,mpireqr2,mpierr)
                Call Mpi_Irecv( id(recvpoint),  nrecv, MPI_PRI,prov,  prov,MpiCube,mpireqr3,mpierr)
            End If
            ! received particles are stored in final arrays at position recvpoint in the arrays.
            ! recvpoint is updated after each reception
            recvpoint=recvpoint+nrecv

            ! wait for send communication to complete and deallocate temporary send arrays
            If(npartvloc(dest+1)/=0) Then
                Call Mpi_Wait(mpireqs1,mpistat,mpierr)
                Deallocate(tmpsendx)
                Call Mpi_Wait(mpireqs2,mpistat,mpierr)
                Deallocate(tmpsendv)
                Call Mpi_Wait(mpireqs3,mpistat,mpierr)
                Deallocate(tmpsendi)
            End If
            ! wait for recv communication to complete
            If(nrecv/=0) Then
                Call Mpi_Wait(mpireqr1,mpistat,mpierr)
                Call Mpi_Wait(mpireqr2,mpistat,mpierr)
                Call Mpi_Wait(mpireqr3,mpistat,mpierr)
            End If
            ! end of loop over the processes
            End Do processus

            ! last step: current process copy the particles he read from files which are located in its subdomain to
            ! the final positions/velocities/id's arrays
            xmin =  CubeCoord(1)      * deltasd
            xmax = (CubeCoord(1) + 1) * deltasd
            ymin =  CubeCoord(2)      * deltasd
            ymax = (CubeCoord(2) + 1) * deltasd
            zmin =  CubeCoord(3)      * deltasd
            zmax = (CubeCoord(3) + 1) * deltasd

            ind = 0
            Do j=1,mynpart
                If(tmpx(1,j)>=xmin .and. tmpx(1,j)<xmax .and. &
                tmpx(2,j)>= ymin .and. tmpx(2,j) < ymax .and. &
                tmpx(3,j)>= zmin .and. tmpx(3,j) < zmax) Then
                x(:,recvpoint+ind) = tmpx(:,j)
                v(:,recvpoint+ind) = tmpv(:,j)
                id(recvpoint+ind)  = tmpi(j)
                ind = ind+1
                End If
            End Do

            ! if current process has not the same number of particles at the end of this subroutine than the number it found
            ! while reading, then there's a bug somewhere
            If(recvpoint+ind /= npartv(procID+1)+1) Then
                Write(tmpstr1,'(I9.9)') recvpoint+ind
                Write(tmpstr2,'(I9.9)') npartv(procID+1)+1
                Call EmergencyStop('Wrong particles number found after send/recv swaps:'//tmpstr1//' ; '//tmpstr2,2)
            End If

            ! nb of particles for current process = number of particles in its subdomain
            mynpart = npartv(procID+1)

            ! deallocate temporary arrays
            Deallocate(tmpx,tmpv,tmpi)
            Deallocate(npartv, npartvloc)

            ! timings if requested
            If(dotimings) Then
            tTailPart = Mpi_Wtime() - timeInt
            tReadRA = Mpi_Wtime() - time0
            End If
            
            
            call mpi_barrier(MPICube,mpierr)
            write(*,*)procid,'mynpart  ',mynpart
            write(*,*)procid,'minmax(x)',minval(x),maxval(x)
            write(*,*)procid,'minmax(v)',minval(v),maxval(v)
            write(*,*)procid,'minmax(i)',minval(id),maxval(id)

    End Subroutine ramses_lecture


    Subroutine emma_lecture()
        Use modvariable
        Use modparam
        Use modmpicom
        Use modtiming
        Implicit none

        ! Local variables
        Character(len=5)               :: ncharcpu
        Character(len=5)               :: ncharsnap
        Character(len=9)               :: tmpstr1, tmpstr2
        Character(len=200)             :: nomfich
        Character(len=4)              :: dumchar
        Character(len=11)              :: grpchar

        Integer(kind=4)                :: i, j, icpu,idim   ! loop variables
        Integer(kind=4)                :: destCoord(3)       ! coords of the destination MPI process in MPI process cart
        Integer(kind=4)                :: nrecv              ! number of elements received in a Mpi_Recv
        Integer(kind=4)                :: recvpoint          ! address of the 1st element received in the local vector
        Integer(kind=4)                :: mynbfile            ! number of RAMSES part files read by local process
        Integer(kind=4)                :: firstp, lastp      ! id of 1st and last RAMSES part file read
        Integer(kind=4), allocatable   :: npartvloc(:), npartv(:)  ! temp and global table of particle numbers for each process
        Integer(kind=4)                :: n_i, n_j, n_k, nsd, ind
        Integer(kind=4)                :: nproc       ! process number  read in RAMSES info file
!!$        Integer(kind=4)                :: ncpu2       ! process number  read in RAMSES part files
!!$        Integer(kind=4)                :: ndim2       ! dimension       read in RAMSES part files
        Integer(kind=4)                :: npartloc    ! particle number read in RAMSES part files
        Integer(kind=4)                :: prov, dest  ! provenance and destination process number for p2p MPI communications
        Integer(kind=4)                :: mpistat(Mpi_Status_Size)   ! status of MPI communication
        Integer(kind=PRI)              :: tmplongint              ! temp integer8 variable
        Integer(kind=PRI), allocatable :: tmpi(:), tmpsendi(:)    ! TYPE VARIABLE EN FONCTION DU NB DE PART
        Integer(kind=4)                :: errcode
!!$        Integer(kind=4)                :: grpnb

!!$        Real(kind=SP), allocatable     :: tmpsimple(:)            ! temporary variable for Ramses v2 output
!!$        Real(kind=DP), allocatable     :: tmpdouble(:)            ! temporary variable for Ramses v3 output
        Real(kind=SP), allocatable     :: tmpsendx(:,:),tmpsendv(:,:)
        Real(kind=SP), allocatable     :: tmpx(:,:), tmpv(:,:)
        Real(kind=SP)                  :: deltasd

        integer(kind=4)                ::nmod

! EMMA addtitional values
        REAL(kind=4)::tloc
        REAL(kind=4)::buffsimple
        Integer(kind=4)::monidx
        Integer(kind=4),allocatable:: npartE(:),npartEloc(:)
        Integer(kind=8),allocatable:: offidx(:)
        
        Write(ncharsnap(1:5),'(I5.5)') snapnum

!END EMMA

        
        ! buffer for Mpi_Pack
        ! 3*bit_size(nres)/8 : 3x Integer(kind=4) nproc, ndim, nres
        h_length = 3* bit_size(nres)/8
        Allocate (header(0:h_length-1))
        h_pos = 0

        ! Timer initialization
        time0 = Mpi_Wtime()

        ! group directory base name: particles files are dispatched in different 'group' directories
        grpchar = 'group_00001'

        ! Process 0 reads parameters in info file
        If(procID == 0) Then
            ! RAMSES v2 particles files
            If(code_index.eq.'RA2') Then
                Print *,'Reading Ramses v2 output...'
                Write(Ulog,*) 'Reading Ramses v2 output...'
                ! RAMSES v3 particles files
            Else if(code_index.eq.'RA3') Then
                Print *,'Reading Ramses v3 output...'
                Write(Ulog,*) 'Reading Ramses v3 output...'
            Else if(code_index.eq.'EMM') Then
                Print *,'Reading EMMA output...'
                Write(Ulog,*) 'Reading EMMA v3 output...'
            End If

            ! info file
            !!!!RY group modif
            !nomfich = trim(pathinput)//'/'//trim(grpchar)//'/'//trim(nameinfo)
!!$            if(grpsize>0) then
!!$               nomfich = trim(pathinput)//'/'//trim(grpchar)//'/'//trim(nameinfo)
!!$            else
!!$               nomfich = trim(pathinput)//'/'//trim(nameinfo)
!!$            endif
            !!!END RY group modif

            nomfich=trim('infosim.txt')

            Print *,'Reading EMMA info file:',trim(nomfich)

            Open(11,file=nomfich, form='formatted', Status='Old', Iostat=errcode)
            If(errcode > 0) Then
                Call EmergencyStop('Error opening '//trim(nameinfo)//' file',errcode)
            End If
            Rewind 11

            ! Read number of particles files, nb of dimensions and level of refinement in the coarse grid (RAMSES is based on AMR)
            Read(11,'(A4,I11)') dumchar,nproc
            Read(11,'(A4,I11)') dumchar,ndim
            Read(11,'(A4,I11)') dumchar,lmin

            Close(11)

            ! nb of grid point in each dimension = 2^lmin
            nres = 2**lmin

            Write(*,*) 'Number of:'
            Write(*,'(A25,I6)') ' - files for each output:',nproc
            Write(*,'(A25,I6)') ' - dimensions:           ',ndim
            Write(*,'(A25,I6)') ' - grid points:          ',nres
            Write(Ulog,*) 'nb_proc = ',nproc,'ndim = ',ndim,'nres = ',nres

            ! Pack these 3 parameters for broadcast
            Call Mpi_Pack(nproc, 1, Mpi_Integer, header, h_length, h_pos, Mpi_Comm_World, mpierr)
            Call Mpi_Pack( ndim, 1, Mpi_Integer, header, h_length, h_pos, Mpi_Comm_World, mpierr)
            Call Mpi_Pack( nres, 1, Mpi_Integer, header, h_length, h_pos, Mpi_Comm_World, mpierr)

        End If

        ! Broadcast
        Call Mpi_Bcast(header,h_length,Mpi_Packed,0,Mpi_Comm_World,mpierr)

        ! Other processes unpack the 3 parameters
        If(procID /= 0) Then

            Call Mpi_Unpack(header, h_length, h_pos, nproc, 1, Mpi_Integer, Mpi_Comm_World, mpierr)
            Call Mpi_Unpack(header, h_length, h_pos,  ndim, 1, Mpi_Integer, Mpi_Comm_World, mpierr)
            Call Mpi_Unpack(header, h_length, h_pos,  nres, 1, Mpi_Integer, Mpi_Comm_World, mpierr)

        End If

        ! total nb of grid points
        ngrid = int(nres,kind=8)**3

        If(allocated(header)) Deallocate(header)

        If(procID==0) Print *,'Reading positions...'
        If(procID==0) Write(Ulog,*) 'Reading positions...'

        ! total nb of particles is initialized to 0
        nptot = 0
        ! each process read mynbfile particles files
        mynbfile = nproc / procNB
        ! index of first and last particles file read by current process
        firstp  = procID * mynbfile + 1
        lastp   = (procID+1) * mynbfile
        
        !!!!When mod(nproc, procNB) .ne. 0!!!RY
        nmod=mod(nproc,procNB)
        if (nmod.ne.0)then
           if(procID.le.nmod-1)then
              mynbfile=mynbfile+1
              firstp  = procID * mynbfile + 1
              lastp   = (procID+1) * mynbfile
           else
              firstp  = procID * mynbfile + 1+nmod
              lastp   = (procID+1) * mynbfile+nmod
           endif
        endif
        !!!!End RY!!!

        ! nb of particles in files read by current process
        mynpart = 0

        ! particles will  be distributed "geographicaly" among the process, i.e. at the end of this subroutine
        ! each process will keep only the particles located in a particular subdomain of the whole simulation domain
        ! allocate array containing the nb of particles that will remain on each process at the end of this subroutine (npartv)
        ! and array containing the nb of particles destined to each process but read by current process (npartvloc)
        Allocate(npartv(procNB))
        Allocate(npartvloc(procNB))


        npartv = 0
        npartvloc = 0


        ! broadcast the number of particles per file
        Allocate(npartEloc(nproc))
        Allocate(npartE(nproc))
        Allocate(offidx(nproc))
        
        npartEloc=0
        npartE=0


        ! nsd = nb of subdomains in each direction, i.e. there are nsd**3 subdomains
        nsd = int(procNB**(1./3.))
        ! deltasd = dimension of a subdomain, i.e. each subdomain is a cube of size deltasd
        deltasd = 1./nsd

        If(procID == 0) Then
            Write(*,*) 'Number of subdomains in each dimension:',nsd
            Write(*,*) 'Size of each subdomain:',deltasd
        End If

        ! definition of the current subdomain by min and max coordinates
        xmin =  CubeCoord(1)      * deltasd
        xmax = (CubeCoord(1) + 1) * deltasd
        ymin =  CubeCoord(2)      * deltasd
        ymax = (CubeCoord(2) + 1) * deltasd
        zmin =  CubeCoord(3)      * deltasd
        zmax = (CubeCoord(3) + 1) * deltasd

        ! loop over the "file/cpu" index characterizing the input particles files to be read by current process
        Do icpu = firstp-1,lastp-1
            Write(ncharcpu(1:5),'(I5.5)') icpu

            !nomfich = trim(pathinput)//'.p'//trim(ncharcpu) ! ORG EMMA
            nomfich = trim(pathinput)//trim(ncharsnap)//'/part_x/x.'//trim(ncharsnap)//'.p'//trim(ncharcpu) ! ALLOCT EMMA

            write(*,*) nomfich

            ! read number of files, dimensions and particles
            open(unit=1,file=nomfich,status='old',access='STREAM',form='UNFORMATTED')
            read(1) npartloc
            read(1) tloc
            close(1)

            ! nb of particles read by current process
            mynpart = mynpart + npartloc

            !store the number of part per file
            npartEloc(icpu+1)=npartloc
        End Do


        If(procID==0) Print *,'npartEloc array read...'


        ! we use a temporary long integer in case of simulations larger than 1024**3
        tmplongint = mynpart
        ! each mynpart is added -the resulst is nptot, total nb of particles - and broadcasted in a single collective communication
        Call Mpi_AllReduce(tmplongint,nptot,1,MPI_PRI,Mpi_Sum,Mpi_Comm_World,mpierr)


        ! reduce npartEloc
        Call Mpi_AllReduce(npartEloc,npartE,nproc,MPI_INTEGER,Mpi_Sum,Mpi_Comm_World,mpierr)

        If(procID==0) Print *,'Reduction done...'

        ! computing index offset

        offidx(1)=0
        do i=2,nproc
           offidx(i)=offidx(i-1)+npartE(i-1)
        enddo

        If(procID==0) Print *,'Offset done...'

        if(offidx(nproc)+npartE(nproc).ne.nptot) then
           write(*,*)'Error in offset idex'
           write(*,*)'Expected',nptot,' obtained',offidx(nproc)+npartE(nproc)
        end if


        If(procID == 0) Then
            Write(* ,*)'There are ',nptot,' DM particles'
            Write(Ulog,*)'There are ',nptot,' DM particles'
        End If


        ! memory allocation for temporary arrays which will contain positions, velocities and id's
        ! of the particles read by current process
        Allocate(tmpx(3,mynpart))
        Allocate(tmpv(3,mynpart))
        Allocate(tmpi(mynpart))

        ! intialization
        tmpx=0.
        mynpart = 0
        monidx=0
        ! Begin timings if requested
        If(dotimings) Then
            Call Mpi_Barrier(MPICube,mpierr)
            timeInt = Mpi_Wtime()
            tInitRead = timeInt - time0
        End If

        ! current process reads its particles files again, but this time it reads position, velocities and id's for each particle
        ! in the file

        Do icpu = firstp-1,lastp-1
           Write(ncharcpu(1:5),'(I5.5)') icpu

           ! READING THE ALLOCT FILES

           ! X
           nomfich = trim(pathinput)//trim(ncharsnap)//'/part_x/x.'//trim(ncharsnap)//'.p'//trim(ncharcpu) ! ALLOCT EMMA
           open(unit=1,file=nomfich,status='old',access='STREAM',form='UNFORMATTED')
           read(1) npartloc
           read(1) tloc

           ! Y
           nomfich = trim(pathinput)//trim(ncharsnap)//'/part_y/y.'//trim(ncharsnap)//'.p'//trim(ncharcpu) ! ALLOCT EMMA
           open(unit=2,file=nomfich,status='old',access='STREAM',form='UNFORMATTED')
           read(2) npartloc
           read(2) tloc

           ! Z
           nomfich = trim(pathinput)//trim(ncharsnap)//'/part_z/z.'//trim(ncharsnap)//'.p'//trim(ncharcpu) ! ALLOCT EMMA
           open(unit=3,file=nomfich,status='old',access='STREAM',form='UNFORMATTED')
           read(3) npartloc
           read(3) tloc


           ! loop over particles
           Do j=1,npartloc
              ! Read positions
              Do idim = 1,ndim
                 Read(idim) buffsimple
                 ! put all positions in tmpx vector
                 if(buffsimple<0.) then
                    buffsimple=buffsimple+1.
                    write(*,*) 'plus'
                 else if(buffsimple>1.) then
                    buffsimple=buffsimple-1.
                    write(*,*) 'minus'
                 end if
                 tmpx(idim,j+mynpart) = buffsimple
              End Do
              
              ! Read velocities in a dummy variable
              !Do idim = 1,ndim
              !   Read(idim+3) buffsimple
              !   ! put all velocities in tmpv vector
              !   tmpv(idim,j+mynpart) = buffsimple
              !End Do
              
              
              ! Read particle id
              !Read(1) buffsimple
              !tmpi(j+mynpart)=buffsimple+1
              tmpi(j+mynpart)=j+offidx(icpu+1)
              
              ! Read dummies
              !Read(1) buffsimple !mass
              !Read(1) buffsimple !epot
              !Read(1) buffsimple !ekin

           enddo

!!$           if(icpu.eq.0) then
!!$              write(*,*) 'Minval--',minval(tmpx(1,:)),minval(tmpx(2,:)),minval(tmpx(3,:))
!!$              write(*,*) 'Maxval--',maxval(tmpx(1,:)),maxval(tmpx(2,:)),maxval(tmpx(3,:))
!!$           end if

           ! loop over the particles read in the current particles file
           Do j = mynpart+1,mynpart+npartloc
              ! particle positions must be >= 0 and <1.0
              ! as the domain is periodic, every position=1.0 is set to 0.0
              Do idim = 1,ndim
                 If(tmpx(idim,j)==1.0) tmpx(idim,j) = 0.0
              End Do
              ! in which subdomain is located this particle?
              n_i = int(tmpx(1,j)/deltasd)
              n_j = int(tmpx(2,j)/deltasd)
              n_k = int(tmpx(3,j)/deltasd)
              ind = nsd**2 *n_i + nsd*n_j + n_k + 1
              ! this means subdomain number 'ind' has one more particle to take care of
              npartvloc(ind) = npartvloc(ind)+1
           End Do
           
           ! close particles file
           Close(1)
           Close(2)
           Close(3)
           ! current process has read mynpart particles
           mynpart = mynpart+npartloc
        End Do

        ! add and broadcast npartvloc so that each process knows how many particles it and its friends will have to analize
        Call Mpi_AllReduce(npartvloc,npartv,procNB,Mpi_Integer,Mpi_Sum,Mpi_Comm_World,mpierr)

        write(*,*) 'ProcID=',procID,sum(npartvloc),sum(npartv),mynpart,monidx,npartvloc(1)

        ! timing if requested
        If(dotimings) Then
            tReadfile = Mpi_Wtime() - timeInt
            timeInt = Mpi_Wtime()
        End If

        ! ------------------------------------------------
        ! sharing out of the particles between processes
        ! ------------------------------------------------

        ! allocate final positions/velocities/id's arrays
        Allocate (x(3,npartv(procID+1)))
        Allocate (v(3,npartv(procID+1)))
        Allocate (id(npartv(procID+1)))

        ! index of the position in the arrays where the data received from other processes must be stored
        recvpoint = 1

        ! loop on processes: each process will received from the other processes the positions/velocities/id's of the particles
        ! it has to analyze
        processus : Do i = 1,procNB - 1
        ! ID of the process current process has to send data to
        dest = mod(procID + i,procNB)
        ! ID of the process current process has to receive data from
        prov = mod(procID + procNB - i, procNB)
        ! Coordinates of dest process in the process grid
        Call Mpi_Cart_coords(MPICube,dest,3,destCoord,mpierr)
        ! current process sends and receives the size of the data it will send/receive
        Call Mpi_Isend(npartvloc(dest+1),1,Mpi_Integer,dest,procID,MpiCube,mpireqs1,mpierr)
        Call Mpi_Irecv(nrecv,1,Mpi_Integer,prov,prov,MpiCube,mpireqr1,mpierr)
        ! current process determines the boundaries of the region it has to send to dest
        xmin =  destCoord(1)      * deltasd
        xmax = (destCoord(1) + 1) * deltasd
        ymin =  destCoord(2)      * deltasd
        ymax = (destCoord(2) + 1) * deltasd
        zmin =  destCoord(3)      * deltasd
        zmax = (destCoord(3) + 1) * deltasd

        ! there are npartvloc(dest+1) particles on the current process that it has to send to dest
        ! if this nb is  greater than 0 then current process proceeds
        If(npartvloc(dest+1)>0) Then
            ! allocation of temporary arrays
            Allocate(tmpsendx(3,npartvloc(dest+1)))
            Allocate(tmpsendv(3,npartvloc(dest+1)))
            Allocate(tmpsendi(  npartvloc(dest+1)))

            ind = 1
            ! loop over local particles
            Do j=1,mynpart
                ! if the particle is located in the cube of dest, then it is added to the send buffers
                If(tmpx(1,j)>=xmin .and. tmpx(1,j)<xmax .and. &
                tmpx(2,j)>= ymin .and. tmpx(2,j) < ymax .and. &
                tmpx(3,j)>= zmin .and. tmpx(3,j) < zmax) Then
                tmpsendx(:,ind) = tmpx(:,j)
                tmpsendv(:,ind) = tmpv(:,j)
                tmpsendi(  ind) = tmpi(j)
                ind=ind+1
                End If
            End Do

            ! if current process has not found the same nb of particles in dest cube as it had found while
            ! reading from files, then there is a bug somewhere
            If(ind/=npartvloc(dest+1)+1) Then
                Print *,'Erreur dans la repartition des particules'
                write(*,*) 'ind=',ind,'npartvloc=',npartvloc(dest+1)+1,xmin,xmax,ymin,ymax,zmin,zmax,dest,procID
                Call Mpi_Finalize(mpierr)
                Stop
            End If

       End If

       ! wait for size of send/recv arrays communication to complete
            Call Mpi_Wait(mpireqs1,mpistat,mpierr)
            Call Mpi_Wait(mpireqr1,mpistat,mpierr)

            ! send particles to dest process if needed
            If(npartvloc(dest+1)/=0) Then
                Call Mpi_Isend(tmpsendx,3*npartvloc(dest+1),Mpi_Real,dest,procID,MpiCube,mpireqs1,mpierr)
                Call Mpi_Isend(tmpsendv,3*npartvloc(dest+1),Mpi_Real,dest,  dest,MpiCube,mpireqs2,mpierr)
                Call Mpi_Isend(tmpsendi,  npartvloc(dest+1), MPI_PRI,dest,procID,MpiCube,mpireqs3,mpierr)
            End If
            ! receive particles from recv process if needed
            If(nrecv/=0) Then
                Call Mpi_Irecv(x(1,recvpoint),3*nrecv,Mpi_Real,prov,  prov,MpiCube,mpireqr1,mpierr)
                Call Mpi_Irecv(v(1,recvpoint),3*nrecv,Mpi_Real,prov,procID,MpiCube,mpireqr2,mpierr)
                Call Mpi_Irecv( id(recvpoint),  nrecv, MPI_PRI,prov,  prov,MpiCube,mpireqr3,mpierr)
            End If
            ! received particles are stored in final arrays at position recvpoint in the arrays.
            ! recvpoint is updated after each reception
            recvpoint=recvpoint+nrecv

            ! wait for send communication to complete and deallocate temporary send arrays
            If(npartvloc(dest+1)/=0) Then
                Call Mpi_Wait(mpireqs1,mpistat,mpierr)
                Deallocate(tmpsendx)
                Call Mpi_Wait(mpireqs2,mpistat,mpierr)
                Deallocate(tmpsendv)
                Call Mpi_Wait(mpireqs3,mpistat,mpierr)
                Deallocate(tmpsendi)
            End If
            ! wait for recv communication to complete
            If(nrecv/=0) Then
                Call Mpi_Wait(mpireqr1,mpistat,mpierr)
                Call Mpi_Wait(mpireqr2,mpistat,mpierr)
                Call Mpi_Wait(mpireqr3,mpistat,mpierr)
            End If
            ! end of loop over the processes
            End Do processus

            ! last step: current process copy the particles he read from files which are located in its subdomain to
            ! the final positions/velocities/id's arrays
            xmin =  CubeCoord(1)      * deltasd
            xmax = (CubeCoord(1) + 1) * deltasd
            ymin =  CubeCoord(2)      * deltasd
            ymax = (CubeCoord(2) + 1) * deltasd
            zmin =  CubeCoord(3)      * deltasd
            zmax = (CubeCoord(3) + 1) * deltasd

            ind = 0
            Do j=1,mynpart
                If(tmpx(1,j)>=xmin .and. tmpx(1,j)<xmax .and. &
                tmpx(2,j)>= ymin .and. tmpx(2,j) < ymax .and. &
                tmpx(3,j)>= zmin .and. tmpx(3,j) < zmax) Then
                x(:,recvpoint+ind) = tmpx(:,j)
                v(:,recvpoint+ind) = tmpv(:,j)
                id(recvpoint+ind)  = tmpi(j)
                ind = ind+1
                End If
            End Do

            ! if current process has not the same number of particles at the end of this subroutine than the number it found
            ! while reading, then there's a bug somewhere
            If(recvpoint+ind /= npartv(procID+1)+1) Then
                Write(tmpstr1,'(I9.9)') recvpoint+ind
                Write(tmpstr2,'(I9.9)') npartv(procID+1)+1
                Call EmergencyStop('Wrong particles number found after send/recv swaps:'//tmpstr1//' ; '//tmpstr2,2)
            End If

            ! nb of particles for current process = number of particles in its subdomain
            mynpart = npartv(procID+1)

            ! deallocate temporary arrays
            Deallocate(tmpx,tmpv,tmpi)
            Deallocate(npartv, npartvloc)

            ! timings if requested
            If(dotimings) Then
            tTailPart = Mpi_Wtime() - timeInt
            tReadRA = Mpi_Wtime() - time0
            End If
            
            
            call mpi_barrier(MPICube,mpierr)
            write(*,*)procid,'mynpart  ',mynpart
            write(*,*)procid,'minmax(x)',minval(x),maxval(x)
            write(*,*)procid,'minmax(v)',minval(v),maxval(v)
            write(*,*)procid,'minmax(i)',minval(id),maxval(id)

    End Subroutine emma_lecture


    !=======================================================================

    ! Write the cube of particles to be analized by current process
    Subroutine outputcube()

        Use modparam
        Use modvariable
        Use modmpicom
        Implicit none

        Integer(kind=PRI) :: i, j
        Character(len=90) :: filecube
        Character(len=5)  :: pid_char
        integer(kind=4)  ::grpnb    !!!!RY group modif 
        Character(len=11):: grpchar !!!!RY group modif 
        character(len=100)::filecmd !!!!RY group modif 

        !!!!RY group modif
        !Write(pid_char(1:5),'(I5.5)') procID
        !filecube = trim(root)//"_cube_"//pid_char
        grpchar = 'group_00001'
        Write(pid_char(1:5),'(I5.5)') procID
        if(grpsize >0.and.nptot>nptot_grp_threshold) then
           grpnb = (procID)/grpsize + 1
           Write(grpchar(7:11),'(I5.5)') grpnb
           filecmd='mkdir '//trim(grpchar)
!!!!          call system(trim(filecmd))     
           call Mpi_Barrier(Mpi_Comm_World,mpierr)
           filecube = trim(grpchar)//'/'//trim(root)//"_cube_"//pid_char
        else
           filecube = trim(root)//"_cube_"//pid_char
        endif
        !END RY group modif

        Open(Unit=Ucub,file=filecube,Form='Unformatted')

        ! write number of particles
        Write(Ucub) int(mynpart,kind=4)
        ! write ID of the FoF process assigned to this cube
        Write(Ucub) procID
        ! write boundary values of the cube
        Write(Ucub) xmin,xmax,ymin,ymax,zmin,zmax
        ! write positions
        Write(Ucub) ((x(j,i),j=1,3),i=1,mynpart)
        ! write velocities
        Write(Ucub) ((v(j,i),j=1,3),i=1,mynpart)
        ! write id's
        Write(Ucub) (id(i),i=1,mynpart)

        Close(Ucub)

    End Subroutine outputcube


    !=======================================================================

    ! write the number of particles (mass) and position of center of mass
    ! for each halo treated by the process whose mass > Mmin
    Subroutine outputmass(ns,smin,nbs,massamas,cdmamas)

        Use modparam
        Use modvariable
        Use modmpicom
        Implicit none

        Integer(kind=4),                   intent(in) :: ns        ! nb of halos with mass > Mmin
        Integer(kind=PRI),                 intent(in) :: smin
        Integer(kind=4)                               :: nbs       ! total nb of halos treated by the process
        Integer(kind=4),  dimension(nbs),  intent(in) :: massamas
        Real(kind=DP),    dimension(3,nbs),intent(in) :: cdmamas

        Character(len=90) :: fileamas
        Character(len=5)  :: pid_char
        Integer(kind=4) :: i
        
        integer(kind=4)  ::grpnb    !!!!RY group modif 
        Character(len=11):: grpchar !!!!RY group modif 
        character(len=100)::filecmd !!!!RY group modif 

        !!!!RY group modif
        !Write(pid_char(1:5),'(I5.5)') procID
        !fileamas = trim(root)//"_masst_"//pid_char
        grpchar = 'group_00001'
        Write(pid_char(1:5),'(I5.5)') procID
        if(grpsize >0.and.nptot>nptot_grp_threshold) then
           grpnb = (procID)/grpsize + 1
           Write(grpchar(7:11),'(I5.5)') grpnb
           filecmd='mkdir '//trim(grpchar)
!!!!           call system(trim(filecmd))
           call Mpi_Barrier(Mpi_Comm_World,mpierr)
           fileamas = trim(grpchar)//'/'//trim(root)//"_masst_"//pid_char
        else
           fileamas = trim(root)//"_masst_"//pid_char
        endif
        !END RY group modif
        
        Open(Umas, File=trim(fileamas), Status='Unknown', Form='Unformatted')
        ! write number of halos with mass > Mmin
        Write(Umas) int(ns,kind=4)
        ! write mass and position for each halo with mass > Mmin
        Do i=1,nbs
            If(massamas(i) >= Mmin ) Write(Umas) int(i+smin-1,kind=8), massamas(i), real(cdmamas(:,i),kind=SP)
        End Do
        Close(Umas)

    End Subroutine outputmass

    !=======================================================================

    ! write position, velocity and ID of every particles contained in every halo with mass > Mmin treated by current process
    Subroutine outputstruct(np,ns,nbs,idf,xf,vf,massamas)

        Use modparam
        Use modvariable
        Use modmpicom
        Implicit none

        Integer(kind=4),                  intent(in) :: np, ns, nbs
        Integer(kind=PRI),dimension(np),  intent(in) :: idf
        Integer(kind=4),  dimension(nbs), intent(in) :: massamas
        Real   (kind=SP), dimension(3,np),intent(in) :: xf, vf


        Integer(kind=4) :: i
        Integer(kind=4) :: j, k, b
        Character(len=90) :: filestrct
        Character(len=5)  :: pid_char
        integer(kind=4)  ::grpnb    !!!!RY group modif 
        Character(len=11):: grpchar !!!!RY group modif 
        character(len=100)::filecmd !!!!RY group modif 

        !!!!RY group modif
        !Write(pid_char(1:5),'(I5.5)') procID
        !filestrct = trim(root)//"_strct_"//pid_char
        grpchar = 'group_00001'
        Write(pid_char(1:5),'(I5.5)') procID
        if(grpsize >0.and.nptot>nptot_grp_threshold) then
           grpnb = (procID)/grpsize + 1
           Write(grpchar(7:11),'(I5.5)') grpnb
           filecmd='mkdir '//trim(grpchar)
!!!!           call system(trim(filecmd))
           call Mpi_Barrier(Mpi_Comm_World,mpierr)
           filestrct = trim(grpchar)//'/'//trim(root)//"_strct_"//pid_char
        else
           filestrct = trim(root)//"_strct_"//pid_char
        endif
        !END RY group modif
        
        Open(Ustr, File=trim(filestrct), Status='Unknown',Form='Unformatted')

        ! number of halo with mass > Mmin
        Write(Ustr) int(ns,kind=4)

        ! b is the 'adress' of the currently examined halo in the positions/velocities/ID's arrays
        b = 1
        ! loop over every halos
        Do i=1, nbs
            ! if mass > Mmin then write
            If(massamas(i) >= Mmin) Then
                ! write mass
                Write(Ustr) massamas(i)
                ! write position of every particles
                Write(Ustr) ((xf(k,j),k=1,3),j=b,b+massamas(i)-1)
                ! write velocities
                Write(Ustr) ((vf(k,j),k=1,3),j=b,b+massamas(i)-1)
                ! write ID's
                Write(Ustr) (idf(j),j=b,b+massamas(i)-1)
                ! advance 'adress'
                b = b+massamas(i)
            Else
                ! advance 'adress' without writing if mass < Mmin
                b = b+massamas(i)
            End If
        End Do

        Close(Ustr)

    End Subroutine outputstruct

End Module modio
