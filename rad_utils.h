REAL RadSolver(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct RGRID *stencil, int stride, REAL dtnew, REAL aexp, int chemonlyad);
void coarse2fine_radlin(struct CELL *cell, struct Rtype *Wi);
void coarse2fine_rad(struct CELL *cell, struct Rtype *Wi);
void sanity_rad(int level,struct RUNPARAMS *param, struct OCT **firstoct, struct CPUINFO *cpu, REAL aexp);
void clean_new_rad(int level,struct RUNPARAMS *param, struct OCT **firstoct, struct CPUINFO *cpu, REAL aexp);
void set_new_rad(int level,struct RUNPARAMS *param, struct OCT **firstoct, struct CPUINFO *cpu, REAL aexp);
