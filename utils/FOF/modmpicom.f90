Module modmpicom

    Use modconst

    Character, dimension(:),allocatable :: header     ! buffer for pack/unpack
    Integer(kind=4) :: mpierr                         ! mpi error code
    Integer(kind=4) :: mpireqs1,mpireqs2,mpireqs3,mpireqs4,mpireqr1,mpireqr2,mpireqr3,mpireqr4 ! mpi non blocking communications request
    Integer(kind=4) :: procID, procNB                 ! process ID and number of processes
    Integer(kind=4) :: h_length, h_pos                ! buffer length and position in the buffer
    Integer(kind=4) :: dims(3)                        ! dimensions of the grid of processes
    Integer(kind=4) :: MPICube                        ! new communicator: cartesian grid of processes
    Integer(kind=4) :: CubeCoord(3)                   ! coordinates of the process in the grid
    Integer(kind=4) :: voisin(6)                      ! process ID of the neighbours of current process in the grid
    Logical :: periods(3)                             ! periodicity of the grid of processes

    Contains

    ! Abort execution when a bug is detected
    Subroutine EmergencyStop(message,errcode)

        Character(len=*), Intent(in) :: message  ! Error message
        Integer, Intent(in) :: errcode           ! Error code: this code is not meaningful

        Write(*,1000) message,procID

        1000 Format('*** ',A,' on process ',I5.5,'. ***')

        Call Mpi_Abort(MPICube,errcode,mpierr)
        Stop

    End Subroutine EmergencyStop

End Module modmpicom
