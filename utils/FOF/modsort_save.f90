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

End Module modsort
