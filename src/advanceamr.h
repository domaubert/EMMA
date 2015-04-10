REAL Advance_level(int level,REAL *adt, struct CPUINFO *cpu, struct RUNPARAMS *param, struct OCT **firstoct,  struct OCT ** lastoct, struct HGRID *stencil,struct STENGRAV *gstencil, struct RGRID *rstencil,int *ndt, int nsteps,REAL tloc);
REAL Advance_level_RAD(int level,REAL dtmax, REAL *adt, struct CPUINFO *cpu, struct RUNPARAMS *param, struct OCT **firstoct,  struct OCT ** lastoct, struct HGRID *stencil, struct STENGRAV *gstencil, struct RGRID *rstencil, int nsteps, REAL tloc,int nrad);


