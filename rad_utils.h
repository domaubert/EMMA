void RadSolver(int level,struct RUNPARAMS *param, struct OCT *** octList, struct OCT ** firstoct, struct CPUINFO *cpu, struct RGRID *stencil, int stride, REAL dtnew, REAL aexp);
void coarse2fine_radlin(struct CELL *cell, struct Rtype *Wi);
void coarse2fine_rad(struct CELL *cell, struct Rtype *Wi);
void sanity_rad(int level,struct RUNPARAMS *param, struct OCT **octList, struct CPUINFO *cpu, REAL aexp);
void clean_new_rad(struct OCT *** octList, int level,struct RUNPARAMS *param, struct OCT **firstoct, struct CPUINFO *cpu, REAL aexp);
