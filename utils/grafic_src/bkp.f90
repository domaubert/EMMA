
!======================================================================================!
program Grafic2Gadget
  !======================================================================================!
  ! Converts ic_velc* files produced by grafic to a gadget IC file IC.gad
  ! CDM only. 
  !
  ! Usage: Grafic2Gadget input_directory output_directory (nfiles)
  ! 
  ! nfiles is optional. if nfiles >1 the gadget file will be splitted in nfiles subboxes.
  ! 
  ! Notes:
  ! * 'Zeldovich Routines' taken from the GRAFIC package written by
  !   E. Bertschinger. 
  ! * To spare memory only velocities are allocated,transformed to positions, 
  !   then transformed back to velocities
  ! * It takes advantage of the layer structure of Grafic to divide the memory 
  !   consumption by NFILES when the output is splitted
  !
  ! (DA) 18/01/2006 Multiple file Version. Ramses old school format removed.
  ! (DA) 18/10/2005 Original Version, translated from a f77 code 
  ! 
  ! D. Aubert (dominique.aubert@cea.fr)
  !======================================================================================!


  implicit none

  integer:: np1,np2,np3
  real,dimension(:,:,:),allocatable:: velx,vely,velz
  real,dimension(:,:,:),allocatable:: posx,posy,posz
  real,dimension(:),allocatable::ind
  
  real:: dx,x1o,x2o,x3o,astart,omegam,omegav,h0,z0,y0,x0
  integer:: i1,i2,i3,i
  integer:: b1,b2,b3
  real:: vfact,fomega,dladt
  external fomega,dladt

  integer*4,dimension(0:5):: npart_out, nall_out
  real*4,dimension(0:5) ::massarr_out(0:5)
  real*4::a_out, redshift_out
  integer*4,dimension(64-6-12-4-8-2-4*2) ::unused_out
  integer*4 ::flag_sfr_out,flag_feedback_out,flag_cooling_out,numfiles_out
  real*4 ::xLbox_out,omega0_out,omegaL_out,hubble_out
  real*8 ::facco

  integer*4::ifile,nmult=1,nmult2
  character*80::fname,nmult2s
  character*2::dummy2
  character*1::dummy1
  real::ii
  integer::count


  integer*4:: narg
  character*80::input,output
  integer*4,dimension(3)::maxerrx
  integer:: iargc

  real::sumv


  narg = iargc()
  IF(narg .LT. 2)THEN
     write(*,*)'You should type: a.out input output'
     write(*,*)'where directory input should contain GRAFIC files'
     write(*,*)'and directory output will receive the file IC.PM'
     write(*,*)'To Split the files in e.g. 8 subfiles, type:'
     write(*,*)'a.out input output 8 '
     STOP
  END IF


  CALL getarg(1,input)
  CALL getarg(2,output)
  CALL getarg(3,nmult2s)

  if(trim(nmult2s).ne.'') then
     read(nmult2s,*) nmult
  end if


  open(18,file=trim(input)//'/ic_velcx',form='unformatted')
  read(18) np1,np2,np3,dx,x1o,x2o,x3o,astart,omegam,omegav,h0
  close(18)


  !*******************
  ! MULTIPLE FILE CASE     
  !*******************



  write(*,*) '**********************************************'
  write(*,*) 'np1=',np1
  write(*,*) 'np2=',np2
  write(*,*) 'np3=',np3
  write(*,*) 'omegam=',omegam
  write(*,*) 'omegav=',omegav
  write(*,*) 'h0=',h0
  write(*,*) 'dx=',dx
  write(*,*) 'np1*dx=',np1*dx
  write(*,*) 'zstart=',1./astart-1.
  write(*,*) '# of files Grafic2Gadget will output: ',nmult 
  write(*,*) '**********************************************'



  open(18,file=trim(input)//'/ic_velcx',form='unformatted')
  open(19,file=trim(input)//'/ic_velcy',form='unformatted')
  open(20,file=trim(input)//'/ic_velcz',form='unformatted')


  !
  ! ALLOCATION
  !

  allocate(velx(np1,np2,np3/nmult))
  allocate(vely(np1,np2,np3/nmult))
  allocate(velz(np1,np2,np3/nmult))

  allocate(posx(np1,np2,np3/nmult))
  allocate(posy(np1,np2,np3/nmult))
  allocate(posz(np1,np2,np3/nmult))
  allocate(ind(np1*np2*np3))

  write(*,*) 'allocation done'

  read(18) np1,np2,np3,dx,x1o,x2o,x3o,astart,omegam,omegav,h0
  read(19) np1,np2,np3,dx,x1o,x2o,x3o,astart,omegam,omegav,h0
  read(20) np1,np2,np3,dx,x1o,x2o,x3o,astart,omegam,omegav,h0

  !
  ! GADGET HEADER SETTINGS
  !

  a_out=astart
  redshift_out=1.d0/a_out-1.d0
  flag_sfr_out=0
  flag_feedback_out=0.
  nall_out=npart_out
  flag_cooling_out=0
  numfiles_out=1
  hubble_out=h0/100.d0
  xLbox_out=1000.*dx*dble(np1)*hubble_out
  omega0_out=omegam
  omegaL_out=omegav
  unused_out=0



  npart_out(:)=0
  npart_out(1)=np1*np2*np3

  massarr_out(:)=0.d0
  !massarr_out(1)=omegam*2.7755d11*(dx*dble(np1)*h0/100.d0)**3/npart_out(1)/1.d10
  massarr_out(1)=1./npart_out(1)
  write(*,*) 'mass per particle = ',massarr_out(1)


  write(*,*) 'Splitting the files.......'

  npart_out(1)=npart_out(1)/nmult
  numfiles_out=nmult


  !
  ! READING VELOCITY FILES + CONVERSION TO POSITIONS
  !
  
  do i3=1,np3/nmult
     read(18) ((velx(i1,i2,i3),i1=1,np1),i2=1,np2)
     read(19) ((vely(i1,i2,i3),i1=1,np1),i2=1,np2)
     read(20) ((velz(i1,i2,i3),i1=1,np1),i2=1,np2)
  end do
  
  write(*,*) 'Generate Zeldovich displacements...'
  
  vfact=fomega(astart,omegam,omegav)*h0*dladt(astart,omegam,omegav)/astart
  
  do i3=1,np3/nmult
     z0=(i3-1)*dx!+(ifile-1)*dx*np3/nmult
     do i2=1,np2
        y0=(i2-1)*dx
        do i1=1,np1
           
           x0=(i1-1)*dx
           
           posx(i1,i2,i3)=x0+velx(i1,i2,i3)/vfact
           posy(i1,i2,i3)=y0+vely(i1,i2,i3)/vfact
           posz(i1,i2,i3)=z0+velz(i1,i2,i3)/vfact
           
           
        end do
     end do
  end do
  
  posx=posx/dx !
  posy=posy/dx ! Length units in DX
  posz=posz/dx!
  
  ! PERIODIC BC
  do i3=1,np3/nmult
     do i2=1,np2
        do i1=1,np1
           if(posx(i1,i2,i3).ge.np1) then
              posx(i1,i2,i3)=posx(i1,i2,i3)-np1
           elseif(posx(i1,i2,i3).lt.0) then
              posx(i1,i2,i3)=posx(i1,i2,i3)+np1
           endif
           
           if(posy(i1,i2,i3).ge.np2) then
              posy(i1,i2,i3)=posy(i1,i2,i3)-np2
           elseif(posy(i1,i2,i3).lt.0) then
              posy(i1,i2,i3)=posy(i1,i2,i3)+np2
           endif
           
           if(posz(i1,i2,i3).ge.np3) then
              posz(i1,i2,i3)=posz(i1,i2,i3)-np3
           elseif(posz(i1,i2,i3).lt.0) then
              posz(i1,i2,i3)=posz(i1,i2,i3)+np3
           endif
           
        end do
     end do
  end do
  
  

  
  !
  ! SWITCH VELOCITIES TO COMOVING MOMENTUM
  !
  
  write(*,*) 'Switch to velocities'
  ii=0.
  count=1
  sumv=0
  do i3=1,np3/nmult
     do i2=1,np2
        do i1=1,np1
           velx(i1,i2,i3)=velx(i1,i2,i3)*astart**1/(dx*h0)
           vely(i1,i2,i3)=vely(i1,i2,i3)*astart**1/(dx*h0)
           velz(i1,i2,i3)=velz(i1,i2,i3)*astart**1/(dx*h0)
           sumv=sumv+sqrt(velx(i1,i2,i3)*velx(i1,i2,i3)+vely(i1,i2,i3)*vely(i1,i2,i3)+velz(i1,i2,i3)*velz(i1,i2,i3))
           ind(count)=ii
           ii=ii+1.
           count=count+1
        end do
     end do
  end do
  
  do ifile=1,nmult
     
     ! GENERATING NAMES
     ! + OPEN FILES
     
     if(ifile-1.ge.10) then
        write(dummy2,'(i2)') ifile-1
        fname=trim(output)//'IC.PM.'//dummy2
     else
        write(dummy1,'(i1)') ifile-1
        fname=trim(output)//'IC.PM.'//dummy1
     end if
     
     write(*,*) 'opening '//trim(fname)
     
     open(100+ifile,file=trim(fname),form='unformatted',status='unknown')
     write(*,*) 'dx*h0=',dx*h0
     write(*,*) 'average vel=',(sumv)/(np1*np2*np3)/astart*(dx*h0)
     write(*,*) 'average vel=',(sumv)/(np1*np2*np3)/astart*(dx*h0)/sqrt(astart)
     
     !
     ! WRITING THE HEADER
     !
     
     write(100+ifile) npart_out(1),massarr_out(1),a_out,xLbox_out,omega0_out,omegaL_out,hubble_out
     
     
     !
     ! WRITING POSITIONS IN DX UNITS
     !
     
     
!     write(100+ifile) (((posx(i1,i2,i3),posy(i1,i2,i3)&
 !         ,posz(i1,i2,i3),i1=1,np1),i2=1,np2),i3=1,np3/nmult)

     write(100+ifile) (((posx(i1,i2,i3),i1=1,np1),i2=1,np2),i3=1,np3/nmult)
     write(100+ifile) (((posy(i1,i2,i3),i1=1,np1),i2=1,np2),i3=1,np3/nmult)
     write(100+ifile) (((posz(i1,i2,i3),i1=1,np1),i2=1,np2),i3=1,np3/nmult)
     ! WRITING VELOCITIES

     write(100+ifile) (((velx(i1,i2,i3),i1=1,np1),i2=1,np2),i3=1,np3/nmult)
     write(100+ifile) (((vely(i1,i2,i3),i1=1,np1),i2=1,np2),i3=1,np3/nmult)
     write(100+ifile) (((velz(i1,i2,i3),i1=1,np1),i2=1,np2),i3=1,np3/nmult)
     
!     write(100+ifile) (((velx(i1,i2,i3),vely(i1,i2,i3),velz(i1,i2,i3),i1=1,np1),i2=1,np2),i3=1,np3/nmult)
     write(100+ifile) (ind(i3),i3=1,np1*np2*np3)
     
     close(100+ifile)
     
     
  end do

  close(18)
  close(19)
  close(20)
  
  
  
end program Grafic2Gadget


!===============================================================================
function dladt(a,omegam,omegav)
!===============================================================================
!  Evaluate dln(a)/dtau for FLRW cosmology.
!  Omegam := Omega today (a=1) in matter.
!  Omegav := Omega today (a=1) in vacuum energy.
!===============================================================================

  implicit none
  real*4:: dladt,a,omegam,omegav,eta
  
  eta=sqrt(omegam/a+omegav*a*a+1.0-omegam-omegav)
  dladt=a*eta
  return
end function dladt


!===============================================================================
function dplus(a,omegam,omegav)
!===============================================================================
!  Evaluate D+(a) (linear growth factor) for FLRW cosmology.
!  Omegam := Omega today (a=1) in matter.
!  Omegav := Omega today (a=1) in vacuum energy.
!===============================================================================
  implicit none
  real*4 ::dplus,a,omegam,omegav,eta
  real*8 :: om,ov,adp 
  real*8 :: ddplus,rombint
  external ddplus,rombint
  common /omegas/om,ov 

  om=omegam
  ov=omegav
  adp=a
  eta=sqrt(omegam/a+omegav*a*a+1.0-omegam-omegav)
  dplus=eta/a*rombint(ddplus,0.0d0,adp,1.0d-8)
  return
end function dplus

!===============================================================================
function fomega(a,omegam,omegav)
!===============================================================================
!  Evaluate f := dlog[D+]/dlog[a] (logarithmic linear growth rate) for
!  lambda+matter-dominated cosmology.
!  Omega0 := Omega today (a=1) in matter only.  Omega_lambda = 1 - Omega0.
!===============================================================================
  implicit none
  real*4 ::fomega,a,omegam,omegav
  real*4 ::omegak,eta
  real*4 ::dplus
  external dplus


  if (omegam.eq.1.0.and.omegav.eq.0.0) then
     fomega=1.0
     return
  end if
  
  omegak=1.0-omegam-omegav
  eta=sqrt(omegam/a+omegav*a*a+omegak)
  fomega=(2.5/dplus(a,omegam,omegav)-1.5*omegam/a-omegak)/(eta*eta)
  return

end function fomega

!===============================================================================
function ddplus(a)
!===============================================================================
  implicit none
  real*8:: a,eta,ddplus
  real*8::omegam,omegav 
  common /omegas/omegam,omegav 

  if (a.eq.0.0d0) then
     ddplus=0.0d0
     return
  end if
  eta=sqrt(omegam/a+omegav*a*a+1.0d0-omegam-omegav)
  ddplus=2.5d0/(eta*eta*eta)
  return
end function ddplus


!===============================================================================
function rombint(f,a,b,tol)
!===============================================================================

  !
  !
  !  Rombint returns the integral from a to b of f(x)dx using Romberg integration.
  !  The method converges provided that f(x) is continuous in (a,b).  The function
  !  f must be double precision and must be declared external in the calling
  !  routine.  tol indicates the desired relative accuracy in the integral.
  !
  
  implicit none
  integer,parameter :: MAXITER=16,MAXJ=5
  real*8 ::a,b,tol,h,g0,fourj,gmax,error,g1,rombint
  real*8, dimension(MAXJ+1)::g
  integer :: nint
  integer :: i,j,jmax,k
  real*8 ::f
  external f
  

  h=0.5d0*(b-a)
  gmax=h*(f(a)+f(b))
  g(1)=gmax
  nint=1
  error=1.0d20
  i=0
  do
     i=i+1
     if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) exit
     !  Calculate next trapezoidal rule approximation to integral.
     g0=0.0d0
     do k=1,nint
        g0=g0+f(a+(k+k-1)*h)
     enddo
     g0=0.5d0*g(1)+h*g0
     h=0.5d0*h
     nint=nint+nint
     jmax=min(i,MAXJ)
     fourj=1.0d0
     do  j=1,jmax
        !  Use Richardson extrapolation.
        fourj=4.0d0*fourj
        g1=g0+(g0-g(j))/(fourj-1.0d0)
        g(j)=g0
        g0=g1
     enddo
     if (abs(g0).gt.tol) then
        error=1.0d0-gmax/g0
     else
        error=gmax
     end if
     gmax=g0
     g(jmax+1)=g0
  enddo
  rombint=g0
  if (i.gt.MAXITER.and.abs(error).gt.tol) then
     write(*,*) 'Rombint failed to converge; integral, error=',rombint,error
  endif
  return
end function rombint
