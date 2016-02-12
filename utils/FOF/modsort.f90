! Sort subroutines.
! Two different algorithms can be used, quick sort and heap sort.
! Heap sort has proven to be more efficient in our use cases.
! This subroutines sort the id, positions and velocities of particles arrays following increasing structure id.
! Details about these sorting algorithms can be found in the reference book used to write them: 
! Introduction to algorithms (Cormen, Leiserson, Rivest, Stein)
Module modsort

  Use modconst

Contains
  
  Recursive Subroutine trirapide(p,r,tref,tid,tx,tv)
    ! Quick sort as described in Introduction to algorithms (Cormen, Leiserson, Rivest, Stein)
    
    Implicit None

    ! Input parameters
    Integer(kind=4), Intent(in) :: p,r
    ! Input parameters modified in subroutine
    Real   (kind=SP), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid ! structure id and particle id
    ! Local parameters
    Integer(kind=4) :: q

    If(p<r) Then
       
       Call partition(tref,p,r,q,tid,tx,tv)
       Call trirapide(p,q-1,tref,tid,tx,tv)
       Call trirapide(q+1,r,tref,tid,tx,tv)

    End If

  End Subroutine trirapide


  Subroutine partition(tref,p,r,q,tid,tx,tv)

    Implicit None

    ! Input parameters
    Integer(kind=4),Intent(in)  :: p, r  ! First and last index of the tables to sort
    ! Input parameters modified in subroutine
    Real   (kind=SP), Intent(inout),dimension(3,*) :: tx,  tv  ! positions and velocities 
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref,tid ! structure id and particle id
    ! Output parameters
    Integer(kind=4), Intent(out) :: q
    ! Local parameters
    Integer(kind=PRI)     :: tmpi,Sref
    Real(SP),dimension(3) :: tmpr
    Integer(kind=4) :: i, j

    Sref  = tref(r)

    i = p-1
    Do j = p, r-1
       If(tref(j) <= Sref) Then
          i = i+1
          tmpi = tref(i)
          tref(i) = tref(j)
          tref(j) = tmpi
          tmpi = tid(i)
          tid(i) = tid(j)
          tid(j) = tmpi
          tmpr(:) = tx(:,i)
          tx(:,i) = tx(:,j)
          tx(:,j) = tmpr(:)
          tmpr(:) = tv(:,i)
          tv(:,i) = tv(:,j)
          tv(:,j) = tmpr(:)
       End If
    End Do

    tmpi = tref(i+1)
    tref(i+1) = tref(r)
    tref(r) = tmpi
    tmpi = tid(i+1)
    tid(i+1) = tid(r)
    tid(r) = tmpi
    tmpr(:) = tx(:,i+1)
    tx(:,i+1) = tx(:,r)
    tx(:,r) = tmpr(:)
    tmpr(:) = tv(:,i+1)
    tv(:,i+1) = tv(:,r)
    tv(:,r) = tmpr(:)

    q = i+1

  End Subroutine partition


  Subroutine tritas(tref,tid,tx,tv,n)

    Use modmpicom
    Implicit None

    ! Input/Output variables
    Real   (kind=SP), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid ! structure id and particle id

    ! Input variable
    Integer(kind=4),Intent(in) :: n

    ! Local variables
    Integer(kind=4) :: i,size

    Call construction(tref,tid,tx,tv,size,n)

    Call Mpi_Barrier(Mpi_Comm_World,mpierr)

!!$    Print *,'PARAMETRE TRI',size,n
!!$    Print *,'Construction du tas terminee'

    Do i = n, 2, -1
       Call echanger(tref,tid,tx,tv,1,i)
       size = size - 1
       Call entasser(tref,tid,tx,tv,size,1)
    End Do

  End Subroutine tritas

  
  Subroutine construction(tref,tid,tx,tv,size,n)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=SP), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid ! structure id and particle id

    !Output variable
    Integer(kind=4),Intent(out)                :: size

    ! Input variable
    Integer(kind=4),Intent(in) :: n

    ! Local variables
    Integer(kind=4) :: i

    size = n
    Do i = n/2, 1, -1
       Call entasser(tref,tid,tx,tv,size,i)
    End Do

  End Subroutine construction


  Recursive Subroutine entasser(tref,tid,tx,tv,size,i)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=SP), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid ! structure id and particle id

    ! Input variable
    Integer(kind=4),Intent(in) :: size,i

    ! Local variables
    Integer(kind=4) :: l, r, max


    l = 2*i
    r = 2*i+1
    
    If( (l <= size) .and. tref(l) > tref(i) ) Then
       max = l
    Else
       max = i
    End If
    
    If( (r <= size) .and. tref(r) > tref(max) ) Then
       max = r
    End If
    
    If( max /= i ) Then
       Call echanger(tref,tid,tx,tv,i,max)
       Call entasser(tref,tid,tx,tv,size,max)
    End If

  End Subroutine entasser


  Subroutine echanger(tref,tid,tx,tv,i,j)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=SP), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid ! structure id and particle id

    ! Input variable
    Integer(kind=4),Intent(in) :: i,j

    ! Local variables
    Integer(kind=PRI) :: tmpi
    Real(kind=SP ), dimension(3) :: tmpr

    tmpi = tref(i)
    tref(i) = tref(j)
    tref(j) = tmpi

    tmpi = tid(i)
    tid(i) = tid(j)
    tid(j) = tmpi

    tmpr(:) = tx(:,i)
    tx(:,i) = tx(:,j)
    tx(:,j) = tmpr(:)

    tmpr(:) = tv(:,i)
    tv(:,i) = tv(:,j)
    tv(:,j) = tmpr(:)


  End Subroutine echanger


End Module modsort
