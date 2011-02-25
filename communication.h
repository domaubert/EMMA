
void  setup_mpi(struct CPUINFO *cpu, struct OCT **firstoct, int levelmax, int levelcoarse, int ngridmax);
void gather_ex(struct CPUINFO *cpu, struct PACKET **sendbuffer, int field);
int gather_ex_part(struct CPUINFO *cpu, struct PART_MPI **psendbuffer,struct PART **lastp);
void gather_mpi(struct CPUINFO *cpu, struct PACKET **sendbuffer, int field);
void scatter_mpi(struct CPUINFO *cpu, struct PACKET **recvbuffer,  int field);
int scatter_mpi_part(struct CPUINFO *cpu, struct PART_MPI **precvbuffer, struct PART **lastp);
void compute_bndkeys(struct CPUINFO *cpu, struct PACKET **recvbuffer);
void  clean_mpibuff(struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer);
void  clean_mpibuff_part(struct CPUINFO *cpu, struct PART_MPI **psendbuffer, struct PART_MPI **precvbuffer);

#ifdef WMPI
void mpi_exchange(struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer, int field);
void mpi_cic_correct(struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer, int field);
int mpi_exchange_part(struct CPUINFO *cpu, struct PART_MPI **psendbuffer, struct PART_MPI **precvbuffer, struct PART **lastpart);
#endif

