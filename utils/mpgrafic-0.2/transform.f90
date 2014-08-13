module transform

  use grafic_types
  IMPLICIT NONE
#ifdef ADD1US
#define  rfftw3d_f77_create_plan  rfftw3d_f77_create_plan_
#define  rfftwnd_f77_destroy_plan rfftwnd_f77_destroy_plan_
#define  rfftwnd_f77_one_real_to_complex rfftwnd_f77_one_real_to_complex_
#define  rfftw3d_f77_mpi_create_plan  rfftw3d_f77_mpi_create_plan_
#define  rfftwnd_f77_mpi_destroy_plan rfftwnd_f77_mpi_destroy_plan_
#define  rfftwnd_f77_mpi rfftwnd_f77_mpi_
#define  rfftwnd_f77_mpi_local_sizes rfftwnd_f77_mpi_local_sizes_
#endif

  interface fft_mpi
     module procedure fft_mpi_single, fft_mpi_double
  end interface

contains

  subroutine fft_forward(l,m,n,input,output)

    integer :: l,m,n
    real*8, dimension(l,m,n), intent(in) :: input
    complex*16, dimension(l/2+1,m,n), intent(out) :: output
    integer*8 :: plan

    call rfftw3d_f77_create_plan(plan,l,m,n,FFTW_REAL_TO_COMPLEX, &
         & FFTW_ESTIMATE)
    call rfftwnd_f77_one_real_to_complex(plan,input,output)
    call rfftwnd_f77_destroy_plan(plan)

  end subroutine fft_forward

  ! This covers both forward and backward ffts, 
  ! depending on plan

  subroutine fft_mpi_double(plan,input,total_local_size)

    ! Arguments
    integer(i8b) :: plan
    integer :: total_local_size
    real(dp), dimension(total_local_size) :: input

    ! Local variables
    integer :: ierr,myid,status
    real*4 :: tstart, tstop
    real*8, allocatable, dimension(:) :: work

    !allocate(work(total_local_size),stat=status)
    !if (status /= 0) then
    !   print*,'Could not allocate work for proc #',myid
    !   stop
    !endif
    allocate(work(1))
!!$    call rfftwnd_f77_mpi(plan,1,input,work, &
!!$         & 0,FFTW_TRANSPOSED_ORDER)
    call rfftwnd_f77_mpi(plan,1,input,work, &
         & 0,FFTW_NORMAL_ORDER) !0 to ignore work, 1 to use it
    deallocate(work)

  end subroutine fft_mpi_double

  subroutine fft_mpi_single(plan,input,total_local_size)

    ! Arguments
    integer(i8b) :: plan
    integer :: total_local_size
    real(sp), dimension(total_local_size) :: input

    ! Local variables
    integer :: ierr,myid,status
    real*4 :: tstart, tstop
    real*4, allocatable, dimension(:) :: work

    !allocate(work(total_local_size),stat=status)
    !if (status /= 0) then
    !   print*,'Could not allocate work for proc #',myid
    !   stop
    !endif
    allocate(work(1))
!!$    call rfftwnd_f77_mpi(plan,1,input,work, &
!!$         & 0,FFTW_TRANSPOSED_ORDER)
    call rfftwnd_f77_mpi(plan,1,input,work, &
         & 0,FFTW_NORMAL_ORDER) !0 to ignore work, 1 to use it
    deallocate(work)

  end subroutine fft_mpi_single

 

end module transform
