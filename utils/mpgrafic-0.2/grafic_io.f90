module grafic_io

  use grafic_types
  implicit none
!  integer, parameter :: nblock = 10
  integer, parameter :: RECORD_FLAG=4 !bytes 
  
  type taille
     integer :: nx
     integer :: ny
     integer :: nz
     real(kind=sp) :: dx
     real(kind=sp) :: lx
     real(kind=sp) :: ly
     real(kind=sp) :: lz
  end type taille

  type cosmo
     real(kind=sp) :: astart
     real(kind=sp) :: omegam
     real(kind=sp) :: omegav
     real(kind=sp) :: h0
  end type cosmo

  interface grafic_write
     module procedure grafic_write_single, grafic_write_double
  end interface

  interface grafic_read
     module procedure grafic_read_single, grafic_read_double
  end interface

contains

  
  ! This routine write a fftw_mpi slice of a cube 
  ! This slice has interior dimension nz, exterior nx.
  ! Beware that interior dimension is padded: nz <- 2(nz/2+1)
  ! Padded region will not be written.
  subroutine grafic_write_double(buffer,local_nz,local_z_start,ny,nx,filename,padding_in, white_in)
    
    ! Arguments
    implicit none
    include 'mpif.h'
    real(dp), dimension(:), intent(in) :: buffer
    integer, intent(in) :: local_nz, local_z_start, ny, nx
    character(len=128), intent(in) :: filename
    !integer, intent(in) :: filenamelen
    logical, optional :: padding_in
    logical, optional :: white_in

    ! Local variables
    real(sp), allocatable, dimension(:) :: tampon
    integer :: record_size
    integer(RECORD_FLAG) :: taille_tampon
    integer :: n2x
    integer(i8b) :: index, index2, offset, length, toto
    integer :: idummy=0
    integer :: i,j,k,myid,ierr
    logical :: padding, white


    padding=.false.
    white=.false.
    if (present(padding_in)) padding = padding_in
    if (present(padding_in)) write(*,*) 'padding is on'
    if (present(white_in)) white = white_in
    ! 2*4 for leading and trailing record sizes 
    n2x = 2*(nx/2+1)
    if (padding) then
       write(*,*) 'padding is true in grafic_write'
       allocate(tampon(ny*n2x))
       taille_tampon = ny*n2x*sp ! in bytes
    else
       allocate(tampon(ny*nx))
       taille_tampon = ny*nx*sp ! in bytes
    endif

    if (white) then
       offset = 4*4 + 2*RECORD_FLAG
    else
       offset = 3*4+8*sp + 2*RECORD_FLAG ! First record
    endif

    if (padding) then
       offset = offset + local_z_start*int(ny*n2x*sp + 2*RECORD_FLAG,8)
    else
       offset = offset + local_z_start*int(ny*nx*sp + 2*RECORD_FLAG,8)
    endif

    do k=1,local_nz
       
       index2=1
       ! This loop to (eventually) get rid of fftw padding zone
       if (padding) then
       do j=1,ny
          do i=1,n2x
             index = ((k-1)*ny+j-1)*n2x+i
             tampon(index2) = real(buffer(index),kind=sp)
             index2 = index2 + 1
          enddo
       enddo
       else
          do j=1,ny
             do i=1,nx
                index = ((k-1)*ny+j-1)*n2x+i
                tampon(index2) = real(buffer(index),kind=sp)
                index2 = index2 + 1
             enddo
          enddo
       endif
       ! First write f... Fortran record length
       length=RECORD_FLAG
       call f77_parallel_write(trim(filename),len_trim(filename), length &
            & ,offset,taille_tampon)
       offset = offset+RECORD_FLAG
       if (padding) then
          length = ny*n2x*sp
       else
          length = ny*nx*sp
       endif
       call f77_parallel_write(trim(filename),len_trim(filename), length, &
            & offset, tampon)
       offset = offset+length
       ! Again ..
       length=RECORD_FLAG
       call f77_parallel_write(trim(filename),len_trim(filename), length &
            & ,offset,taille_tampon)
       offset = offset+RECORD_FLAG

    enddo
    deallocate(tampon)

  end subroutine grafic_write_double
       
  subroutine grafic_write_single(buffer,local_nz,local_z_start,ny,nx,filename,padding_in, white_in)
    
    ! Arguments
    implicit none
    include 'mpif.h'
    real(sp), dimension(:), intent(in) :: buffer
    integer, intent(in) :: local_nz, local_z_start, ny, nx
    character(len=128), intent(in) :: filename
    !integer, intent(in) :: filenamelen
    logical, optional :: padding_in
    logical, optional :: white_in


    ! Local variables
    real(sp), allocatable, dimension(:) :: tampon
    integer :: record_size
    integer(RECORD_FLAG) :: taille_tampon
    integer :: n2x
    integer(i8b) :: index, index2, offset, length, toto
    integer :: idummy=0
    integer :: i,j,k,myid,ierr
    logical :: padding, white


    padding=.false.
    white=.false.
    if (present(padding_in)) padding = padding_in
    if (present(padding_in)) write(*,*) 'padding is on'
    if (present(white_in)) white = white_in
    ! 2*4 for leading and trailing record sizes 
    n2x = 2*(nx/2+1)
    if (padding) then
       write(*,*) 'padding is true in grafic_write'
       allocate(tampon(ny*n2x))
       taille_tampon = ny*n2x*sp ! in bytes
    else
       allocate(tampon(ny*nx))
       taille_tampon = ny*nx*sp ! in bytes
    endif

    if (white) then
       offset = 4*4 + 2*RECORD_FLAG
    else
       offset = 3*4+8*sp + 2*RECORD_FLAG ! First record
    endif

    if (padding) then
       offset = offset + local_z_start*int(ny*n2x*sp + 2*RECORD_FLAG,8)
    else
       offset = offset + local_z_start*int(ny*nx*sp + 2*RECORD_FLAG,8)
    endif

    do k=1,local_nz
       
       index2=1
       ! This loop to (eventually) get rid of fftw padding zone
       if (padding) then
       do j=1,ny
          do i=1,n2x
             index = ((k-1)*ny+j-1)*n2x+i
             tampon(index2) = real(buffer(index),kind=sp)
             index2 = index2 + 1
          enddo
       enddo
       else
          do j=1,ny
             do i=1,nx
                index = ((k-1)*ny+j-1)*n2x+i
                tampon(index2) = real(buffer(index),kind=sp)
                index2 = index2 + 1
             enddo
          enddo
       endif
       ! First write f... Fortran record length
       length=RECORD_FLAG
       call f77_parallel_write(trim(filename),len_trim(filename), length &
            & ,offset,taille_tampon)
       offset = offset+RECORD_FLAG
       if (padding) then
          length = ny*n2x*sp
       else
          length = ny*nx*sp
       endif
       call f77_parallel_write(trim(filename),len_trim(filename), length, &
            & offset, tampon)
       offset = offset+length
       ! Again ..
       length=RECORD_FLAG
       call f77_parallel_write(trim(filename),len_trim(filename), length &
            & ,offset,taille_tampon)
       offset = offset+RECORD_FLAG

    enddo
    deallocate(tampon)

  end subroutine grafic_write_single

  subroutine grafic_read_double(buffer,local_nz,local_z_start,ny,nx,filename, padding_in, white_in)

    ! Arguments
    implicit none
    include 'mpif.h'
    integer, intent(in) :: local_nz, local_z_start, ny, nx
    real(dp), dimension(local_nz*ny*2*(nx/2+1)), intent(out) :: buffer
    character(len=128), intent(in) :: filename
    logical, optional :: padding_in ! Read padded zone or not
    logical, optional :: white_in ! Different header for white noise file
    ! Local variables
    real(sp), allocatable, dimension(:) :: tampon
    integer :: record_size
    integer(RECORD_FLAG) :: taille_tampon
    integer :: n2x
    integer(i8b) :: index, index2, offset, length
    integer :: idummy=0
    integer :: i,j,k
    integer :: myid,ierr
    logical :: padding
    logical :: white

    padding = .false.
    white = .false.
    if (present(padding_in)) padding = padding_in
    if (present(white_in)) white = white_in
    ! 2*4 for leading and trailing record sizes 
    n2x = 2*(nx/2+1)
    if (padding) then
       allocate(tampon(ny*n2x))
    else 
       allocate(tampon(ny*nx))
    endif
    taille_tampon = 0

    if (white) then
       offset = 4*4 + 2*RECORD_FLAG
    else
       offset = 3*4+8*sp + 2*RECORD_FLAG ! First record
    endif
    if (padding) then
       offset = offset + local_z_start*int(ny*n2x*sp + 2*RECORD_FLAG,8)
    else 
       offset = offset + local_z_start*int(ny*nx*sp+2*RECORD_FLAG,8)
    endif

    do k=1,local_nz
       
       length = RECORD_FLAG
       ! First read f... Fortran record length
       call f77_parallel_read(trim(filename),len_trim(filename), length &
            & ,offset,taille_tampon)
       if (padding) then
           if (taille_tampon /= ny*n2x*sp) then
             print*,'Record size ',taille_tampon,' is different from expected', &
                  & ny*n2x*sp
             stop
          endif
       else
           if (taille_tampon /= ny*nx*sp) then
             print*,'Record size ',taille_tampon,' is different from expected', &
                  & ny*nx*sp
             stop
          endif
       endif
       offset = offset+RECORD_FLAG

       if (padding) then
          length = ny*n2x*sp
       else 
          length = ny*nx*sp
       endif
       call f77_parallel_read(trim(filename),len_trim(filename), length, &
            & offset, tampon)
       offset = offset+length
       ! Again ..
       length=RECORD_FLAG
       call f77_parallel_read(trim(filename),len_trim(filename), length &
            & ,offset,taille_tampon)
       offset = offset+RECORD_FLAG

       ! This loop to (eventually) get rid of fftw padding zone
       index2=1
       if (padding) then
          do j=1,ny
             do i=1,n2x ! Read padding (Nyquist) : reading in k space ...
                index=((k-1)*ny+j-1)*n2x+i
                buffer(index) = real(tampon(index2),kind=dp)
                index2 = index2 + 1
             enddo
          enddo
       else
          do j=1,ny
             do i=1,nx ! Do not read padding : reading in real space ...
                index = ((k-1)*ny+j-1)*n2x+i
                buffer(index) = real(tampon(index2),kind=dp)
                index2 = index2 + 1
             enddo
          enddo
       endif

    enddo
    deallocate(tampon)

  end subroutine grafic_read_double

  subroutine grafic_read_single(buffer,local_nz,local_z_start,ny,nx,filename, padding_in, white_in)

    ! Arguments
    implicit none
    include 'mpif.h'
    integer, intent(in) :: local_nz, local_z_start, ny, nx
    real(sp), dimension(local_nz*ny*2*(nx/2+1)), intent(out) :: buffer
    character(len=128), intent(in) :: filename
    logical, optional :: padding_in ! Read padded zone or not
    logical, optional :: white_in ! Different header for white noise file
    ! Local variables
    real(sp), allocatable, dimension(:) :: tampon
    integer :: record_size
    integer(RECORD_FLAG) :: taille_tampon
    integer :: n2x
    integer(i8b) :: index, index2, offset, length
    integer :: idummy=0
    integer :: i,j,k
    integer :: myid,ierr
    logical :: padding
    logical :: white

    padding = .false.
    white = .false.
    if (present(padding_in)) padding = padding_in
    if (present(white_in)) white = white_in
    ! 2*4 for leading and trailing record sizes 
    n2x = 2*(nx/2+1)
    if (padding) then
       allocate(tampon(ny*n2x))
    else 
       allocate(tampon(ny*nx))
    endif
    taille_tampon = 0

    if (white) then
       offset = 4*4 + 2*RECORD_FLAG
    else
       offset = 3*4+8*sp + 2*RECORD_FLAG ! First record
    endif
    if (padding) then
       offset = offset + local_z_start*int(ny*n2x*sp + 2*RECORD_FLAG,8)
    else 
       offset = offset + local_z_start*int(ny*nx*sp+2*RECORD_FLAG,8)
    endif

    do k=1,local_nz
       
       length = RECORD_FLAG
       ! First read f... Fortran record length
       call f77_parallel_read(trim(filename),len_trim(filename), length &
            & ,offset,taille_tampon)
       if (padding) then
           if (taille_tampon /= ny*n2x*sp) then
             print*,'Record size ',taille_tampon,' is different from expected', &
                  & ny*n2x*sp
             stop
          endif
       else
           if (taille_tampon /= ny*nx*sp) then
             print*,'Record size ',taille_tampon,' is different from expected', &
                  & ny*nx*sp
             stop
          endif
       endif
       offset = offset+RECORD_FLAG

       if (padding) then
          length = ny*n2x*sp
       else 
          length = ny*nx*sp
       endif
       call f77_parallel_read(trim(filename),len_trim(filename), length, &
            & offset, tampon)
       offset = offset+length
       ! Again ..
       length=RECORD_FLAG
       call f77_parallel_read(trim(filename),len_trim(filename), length &
            & ,offset,taille_tampon)
       offset = offset+RECORD_FLAG

       ! This loop to (eventually) get rid of fftw padding zone
       index2=1
       if (padding) then
          do j=1,ny
             do i=1,n2x ! Read padding (Nyquist) : reading in k space ...
                index=((k-1)*ny+j-1)*n2x+i
                buffer(index) = real(tampon(index2),kind=sp)
                index2 = index2 + 1
             enddo
          enddo
       else
          do j=1,ny
             do i=1,nx ! Do not read padding : reading in real space ...
                index = ((k-1)*ny+j-1)*n2x+i
                buffer(index) = real(tampon(index2),kind=sp)
                index2 = index2 + 1
             enddo
          enddo
       endif

    enddo
    deallocate(tampon)

  end subroutine grafic_read_single




  subroutine grafic_read_header(filename,head_taille,head_cosmo)
    
    type(taille), intent(out) :: head_taille
    type(cosmo), intent(out) :: head_cosmo
    character(len=128), intent(in) :: filename
    logical :: ok

    inquire(file=filename,exist=ok)
    if (.not. ok) then
       print*,'File ',trim(filename),' does not exist, aborting.'
       stop
    endif

    open(2,file=filename,form='unformatted',status='old')
    read(2) head_taille%nx, head_taille%ny, head_taille%nz, head_taille%dx, &
         head_taille%lx, head_taille%ly, head_taille%lz, head_cosmo%astart, &
         head_cosmo%omegam, head_cosmo%omegav, head_cosmo%h0
    close(2)

  end subroutine grafic_read_header

  subroutine grafic_read_header_white(filename,nx,ny,nz,iseed)
    
    integer, intent(out) :: nx,ny,nz,iseed
    character(len=128), intent(in) :: filename
    logical :: ok

    inquire(file=filename,exist=ok)
    if (.not. ok) then
       print*,'File ',trim(filename),' does not exist, aborting.'
       stop
    endif

    open(2,file=filename,form='unformatted',status='old')
    read(2) nx, ny, nz, iseed
    close(2)

  end subroutine grafic_read_header_white
  
  subroutine grafic_write_header(filename,head_taille,head_cosmo)
    
    type(taille), intent(in) :: head_taille
    type(cosmo), intent(in) :: head_cosmo
    character(len=128), intent(in) :: filename

    open(2,file=filename,form='unformatted',status='unknown')
    write(2) head_taille%nx, head_taille%ny, head_taille%nz, head_taille%dx, &
         head_taille%lx, head_taille%ly, head_taille%lz, head_cosmo%astart, &
         head_cosmo%omegam, head_cosmo%omegav, head_cosmo%h0
    close(2)

  end subroutine grafic_write_header

  subroutine grafic_write_header_white(filename,nx,ny,nz,iseed)
    
    integer, intent(in) :: nx,ny,nz,iseed
    character(len=128), intent(in) :: filename

    open(2,file=filename,form='unformatted',status='unknown')
    write(2) nx, ny, nz, iseed
    close(2)

  end subroutine grafic_write_header_white

       
end module grafic_io


  
    
