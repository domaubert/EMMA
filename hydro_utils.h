#ifdef WHYDRO
int hydrosolve(struct MULTIVECT *data, int level, int curcpu, int nread,int stride,REAL dx, REAL dt);
#endif


int hydroS(struct HGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, REAL dt);
int hydroM(struct HGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, REAL dt);
REAL comptstep_hydro(int levelcoarse,int levelmax,struct OCT** firstoct, REAL fa, REAL fa2, struct CPUINFO* cpu, REAL tmax);
REAL comptstep_ff(int levelcoarse,int levelmax,struct OCT** firstoct, REAL aexp, struct CPUINFO* cpu, REAL tmax);
void correct_grav_hydro(struct OCT *octstart, struct CPUINFO *cpu, REAL dt);
REAL comptstep_force(int levelcoarse,int levelmax,struct OCT** firstoct, REAL aexp, struct CPUINFO* cpu, REAL tmax);
void grav_correction(int level,struct RUNPARAMS *param, struct OCT ** firstoct, struct CPUINFO *cpu, REAL dt);
void hydro(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct HGRID *stencil, int stride, REAL dtnew);
REAL L_comptstep_ff(int level,struct RUNPARAMS *param,struct OCT** firstoct, REAL aexp, struct CPUINFO* cpu, REAL tmax);
REAL L_comptstep_hydro(int level, struct RUNPARAMS *param,struct OCT** firstoct, REAL fa, REAL fa2, struct CPUINFO* cpu, REAL tmax);
void HydroSolver(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct HGRID *stencil, int stride, REAL dtnew);
void clean_new_hydro(int level,struct RUNPARAMS *param, struct OCT **firstoct, struct CPUINFO *cpu);
