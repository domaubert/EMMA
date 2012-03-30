#ifdef WHYDRO
int hydrosolve(struct MULTIVECT *data, int level, int curcpu, int nread,int stride,float dx, float dt);
#endif

#ifdef WHYDRO2
#endif

int hydroS(struct HGRID *stencil, int level, int curcpu, int nread,int stride,float dx, float dt);
float comptstep_hydro(int levelcoarse,int levelmax,struct OCT** firstoct, float fa, float fa2, struct CPUINFO* cpu, float tmax);
float comptstep_ff(int levelcoarse,int levelmax,struct OCT** firstoct, float aexp, struct CPUINFO* cpu, float tmax);
void correct_grav_hydro(struct OCT *octstart, struct CPUINFO *cpu, float dt);
