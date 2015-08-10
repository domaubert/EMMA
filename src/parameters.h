void GetParameters(char *fparam, struct RUNPARAMS *param);
void dumpHeader(struct RUNPARAMS *param, struct CPUINFO *cpu,char *fparam);
void dumpStepInfo(struct OCT **firstoct, struct RUNPARAMS *param, struct CPUINFO *cpu, int nsteps,REAL dt,REAL t);
void copy_param(const char *folder);


