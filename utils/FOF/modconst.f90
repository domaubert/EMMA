Module modconst

  Use mpi

#ifdef LONGREAL
  Integer, parameter :: PR=8
#else
  Integer, parameter :: PR=4
#endif

! Define LONGINT in makefile to analyze simulations with more than 2**31 particles.
#ifdef LONGINT
  Integer, parameter :: PRI = 8
  Integer, parameter :: MPI_PRI = Mpi_Integer8
#else
  Integer, parameter :: PRI = 4
  Integer, parameter :: MPI_PRI = Mpi_Integer
#endif


  Integer, parameter :: DP = kind(1.d0)
  Integer, parameter :: SP = kind(1.e0)

  ! Output Units 
  Integer, parameter :: Ulog=50, Ucub=51, Umas=52, Ustr=53, Uopa=54


End Module modconst
