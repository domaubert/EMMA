void breakmpi();
REAL multicheck(struct OCT **firstoct,int *npart,int , int levelmax, int rank, struct CPUINFO *cpu,struct RUNPARAMS *param,int label);
void myradixsort(int *a,int n);
void grid_census(struct RUNPARAMS *param, struct CPUINFO *cpu);
void checkMtot(struct OCT **firstoct, struct RUNPARAMS *param, struct CPUINFO *cpu);
void myradixsort(int *a,int n);

REAL rdm(REAL a, REAL b);
int gpoiss(REAL lambda);
REAL a2t(struct RUNPARAMS *param, REAL az );