! arrays of positions/velocities/ids
Module modvariable

  Use modconst
  Real(   kind=SP ), dimension(:,:), allocatable :: x   ! particles positions - Always simple precision
  Real(   kind=SP ), dimension(:,:), allocatable :: v   ! particles velocities - Always simple precision
  Integer(kind=PRI), dimension(:),   allocatable :: id  ! particles Ramses Id's - Integer8 if particle number 
                                                        ! exceeds 2^31-1
End Module modvariable
