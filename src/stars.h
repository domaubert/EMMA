int testCond(struct CELL *cell, struct RUNPARAMS *param, REAL aexp, int level);
void conserveField(struct Wtype *field, struct RUNPARAMS *param, struct PART *star, REAL dx, REAL aexp, REAL mstar);
void Stars(struct RUNPARAMS *param, struct CPUINFO *cpu, REAL dt,REAL aexp, int level, int is);
