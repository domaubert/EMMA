void RadSolver(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct RGRID *stencil, int stride, REAL dtnew);
void coarse2fine_radlin(struct CELL *cell, struct Rtype *Wi);
void coarse2fine_rad(struct CELL *cell, struct Rtype *Wi);
