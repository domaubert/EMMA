struct PART* findlastpart(struct PART* phead);
int countpart(struct PART* phead);
void countpartDM(struct CELL* cell, int *npart);
struct PART* modifpospart(struct PART* phead, REAL len, int dir);
REAL movepart(int levelcoarse,int levelmax,struct OCT** firstoct, REAL dt, struct CPUINFO *cpu);
void  partcellreorg(int levelcoarse,int levelmax,struct OCT **firstoct);
void forcevel(int levelcoarse,int levelmax,struct OCT **firstoct, REAL **vcomp,int stride,REAL dt, struct CPUINFO *cpu, struct PACKET **sendbuffer,struct PACKET **recvbuffer);
REAL comptstep(int levelcoarse,int levelmax,struct OCT** firstoct, REAL fa, REAL fa2, struct CPUINFO* cpu,REAL);
REAL L_comptstep(int level,struct RUNPARAMS *param,struct OCT** firstoct, REAL fa, REAL fa2, struct CPUINFO* cpu, REAL tmax);
REAL L_movepart(int level,struct OCT** firstoct, REAL *adt, int is, struct CPUINFO* cpu);
void L_accelpart(int level,struct OCT **firstoct, REAL *adt, int is, struct CPUINFO *cpu);
void L_partcellreorg(int level,struct OCT **firstoct);
void L_levpart(int level,struct OCT** firstoct,int is);
REAL L_egypart(int level,struct OCT **firstoct);
void egypart(struct CPUINFO *cpu, REAL *ekintot, REAL *epottot,struct RUNPARAMS *param, REAL tsim);
void L_reset_is_part(int level,struct OCT** firstoct);

void printPart(struct PART* part);
int checkPartNan(struct PART* part);

