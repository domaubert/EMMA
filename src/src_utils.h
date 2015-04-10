void setUVvalue(struct RUNPARAMS *param, REAL aexp);
int FillRad(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, int init, REAL aexp, REAL tloc);
void homosource(struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, int levext);
void cleansource(struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu);
void setUVBKG(struct RUNPARAMS *param, char *fname);
