! Input parameters
Module modparam

  Use modconst
  Character(len=90) :: root               ! base name for output
  Character(len=90) :: pathinput          ! path to the output groups
  Character(len=90) :: nameinfo           ! info file name
  Character(len=90) :: namepart           ! particles files base name
  Character(len=3)  :: code_index         ! input file type: RA2= RAMSES v2, RA3= RAMSES v3
  Integer(kind=4)   :: snapnum            ! snapnumber
  Integer(kind=4)   :: Mmin, Mmax         ! min and max structure mass
  Integer(kind=4)   :: grpsize            ! number of particles files in each group_N directory
  Real(kind=SP)     :: perco              ! percolation parameter for Friends of Friends halo detection
  Logical           :: outcube            ! should there be an output after the particles reading/sharing?
  Logical           :: star               ! is there information about stars in the RAMSES files?
  Logical           :: metal              ! and about metallicity?
  Logical           :: dofof              ! should the structures be detected?
  Logical           :: readfromcube       ! should the particles be read from cube files?
  Logical           :: dotimings          ! should there be timings?

End Module modparam
