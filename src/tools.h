void breakmpi();
REAL multicheck(struct OCT **firstoct,int *npart,int , int levelmax, int rank, struct CPUINFO *cpu,struct RUNPARAMS *param,int label);
void myradixsort(int *a,int n);
void grid_census(struct RUNPARAMS *param, struct CPUINFO *cpu);
void checkMtot(struct OCT **firstoct, struct RUNPARAMS *param, struct CPUINFO *cpu);
void myradixsort(int *a,int n);

double rdm(double a, double b);
unsigned int gpoiss(double lambda);
