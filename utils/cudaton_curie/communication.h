extern "C" int getlocalrank(void);
extern "C" void switchbuff(float *buff, int neighbor, int ndata);
extern "C" void mpisynch(void);
extern "C" void topo_cartesian(int rank, int *dims, int *coords, int *neigh);
extern "C" void get_elapsed(double *);
extern "C" double mpireducemax(double *x);
