void save_amr(char filename[], struct OCT **firstoct,REAL tsim, REAL tinit,int nsteps, int ndumps, struct RUNPARAMS *param, struct CPUINFO *cpu, struct PART *proot, REAL *adt);
struct OCT * restore_amr(char filename[], struct OCT **firstoct,struct OCT **lastoct, REAL *tsim, REAL *tinit, int *nsteps, int *ndumps,struct RUNPARAMS *param, struct CPUINFO *cpu, struct PART *proot, REAL *adt, struct CELL *root);
void save_part(char filename[],struct OCT **firstoct, int levelcoarse, int levelmax, REAL tsim, struct CPUINFO *cpu, struct PART* proot);
struct PART * restore_part(char filename[], struct OCT **firstoct, REAL *tsim, struct RUNPARAMS *param, struct CPUINFO *cpu, struct PART *proot);

