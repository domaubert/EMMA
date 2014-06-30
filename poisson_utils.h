void coarse2fine_grav(struct CELL *cell, struct Gtype *Wi);
struct OCT *gatherstencilgrav(struct OCT ** octList, int iOct,struct OCT *octstart, struct STENGRAV *stencil, int stride, struct CPUINFO *cpu, int *nread);
struct OCT *scatterstencilgrav(struct OCT ** octList, int iOct,struct OCT *octstart, struct STENGRAV *stencil, int nread,int stride, struct CPUINFO *cpu);
int PoissonJacobi(struct OCT ** octList, int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct STENGRAV *stencil, int stride,REAL tsim);
REAL PoissonMgrid(struct OCT *** octList, int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct STENGRAV *stencil, int stride, REAL tsim);
void PoissonForce(struct OCT ** octList, int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct STENGRAV *stencil, int stride, REAL tsim);
int FillDens(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu);
int PoissonSolver(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct STENGRAV *stencil, int stride, REAL aexp);
