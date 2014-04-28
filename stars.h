void createStars(struct OCT **firstoct, struct RUNPARAMS *param, struct CPUINFO *cpu, REAL dt,REAL aexp, int level, int is);
int addStar(struct CELL * cell, int level, REAL xc, REAL yc, REAL zc, struct CPUINFO *cpu, REAL dt,struct RUNPARAMS *param, REAL aexp, REAL drho, int is, REAL dthydro, int nstars);
int  testCond(struct CELL *curcell,REAL dt, REAL dx, struct RUNPARAMS *param, REAL aexp, int level);
void initStar(struct CELL *cell, struct PART *star,struct RUNPARAMS *param, int level,  REAL m ,REAL xc, REAL yc, REAL zc, int idx, REAL aexp, int is, REAL dthydro);
REAL getdrho(struct CELL *cell, REAL dt, REAL aexp);
void removeMfromgas(struct CELL * cell,struct PART *star, REAL m);
