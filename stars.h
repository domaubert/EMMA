void createStars(struct OCT **firstoct, struct RUNPARAMS *param, struct CPUINFO *cpu, REAL dt,REAL aexp, int level);
int addStar(struct CELL * cell, int level, REAL xc, REAL yc, REAL zc, struct CPUINFO *cpu, REAL dt,struct RUNPARAMS *param, REAL aexp, REAL drho);
int  testCond(struct CELL *curcell,REAL dt, REAL dx, struct RUNPARAMS *param, REAL aexp);
void initStar(struct CELL *cell, struct PART *star, int level,  REAL m ,REAL xc, REAL yc, REAL zc, int idx, REAL aexp);
REAL getdrho(struct CELL *cell, REAL dt, REAL aexp);
void removeMfromgas(struct CELL * cell,struct PART *star, REAL m);
int  countStars(struct CELL * cell);
void accretion(struct CELL * cell, struct RUNPARAMS *param, REAL dt, int nstarsincell, REAL dx, REAL aexp);

