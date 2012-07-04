
void coarse2fine_grav(struct CELL *cell, struct Gtype *Wi);
struct OCT *gatherstencilgrav(struct OCT *octstart, struct HGRID *stencil, int stride, struct CPUINFO *cpu, int *nread);
struct OCT *scatterstencilgrav(struct OCT *octstart, struct HGRID *stencil, int stride, struct CPUINFO *cpu);
int PoissonJacobi(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct HGRID *stencil, int stride,REAL tsim);
REAL PoissonMgrid(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct HGRID *stencil, int stride, REAL tsim);
int PoissonForce(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct HGRID *stencil, int stride);
int FillDens(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu);
