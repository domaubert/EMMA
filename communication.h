
void  setup_mpi(struct CPUINFO *cpu, struct OCT **firstoct, int levelmax, int levelcoarse, int ngridmax, int loadb);
void gather_ex(struct CPUINFO *cpu, struct PACKET **sendbuffer, int field);
void gather_ex_part(struct CPUINFO *cpu, struct PART_MPI **psendbuffer);
void gather_mpi(struct CPUINFO *cpu, struct PACKET **sendbuffer, int field);
void scatter_mpi(struct CPUINFO *cpu, struct PACKET **recvbuffer,  int field);
void scatter_mpi_part(struct CPUINFO *cpu, struct PART_MPI **precvbuffer);
void compute_bndkeys(struct CPUINFO *cpu, struct PACKET **recvbuffer);
void clean_mpibuff(struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer);
void clean_mpibuff_part(struct CPUINFO *cpu, struct PART_MPI **psendbuffer, struct PART_MPI **precvbuffer);

#ifdef WMPI
void mpi_exchange(struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer, int field, int cmp_key);
void mpi_cic_correct(struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer, int field);
void mpi_exchange_part(struct CPUINFO *cpu, struct PART_MPI **psendbuffer, struct PART_MPI **precvbuffer, int *delta,int);
void mpi_exchange_level(struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer, int field, int cmp_keys, int level);
void mpi_exchange_hydro(struct CPUINFO *cpu, struct HYDRO_MPI **sendbuffer, struct HYDRO_MPI **recvbuffer, int cmp_keys);
void mpi_hydro_correct(struct CPUINFO *cpu, struct HYDRO_MPI **sendbuffer, struct HYDRO_MPI **recvbuffer,int level);
void mpi_exchange_rad_level(struct CPUINFO *cpu, struct RAD_MPI **sendbuffer, struct RAD_MPI **recvbuffer, int cmp_keys,int level);
void mpi_exchange_hydro_level(struct CPUINFO *cpu, struct HYDRO_MPI **sendbuffer, struct HYDRO_MPI **recvbuffer, int cmp_keys, int level);
#endif

