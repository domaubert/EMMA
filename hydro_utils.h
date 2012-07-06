#ifdef WHYDRO
int hydrosolve(struct MULTIVECT *data, int level, int curcpu, int nread,int stride,REAL dx, REAL dt);
#endif

#ifdef WHYDRO2
#endif

int hydroS(struct HGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, REAL dt);
int hydroM(struct HGRID *stencil, int level, int curcpu, int nread,int stride,REAL dx, REAL dt);
REAL comptstep_hydro(int levelcoarse,int levelmax,struct OCT** firstoct, REAL fa, REAL fa2, struct CPUINFO* cpu, REAL tmax);
REAL comptstep_ff(int levelcoarse,int levelmax,struct OCT** firstoct, REAL aexp, struct CPUINFO* cpu, REAL tmax);
void correct_grav_hydro(struct OCT *octstart, struct CPUINFO *cpu, REAL dt);
REAL comptstep_force(int levelcoarse,int levelmax,struct OCT** firstoct, REAL aexp, struct CPUINFO* cpu, REAL tmax);
