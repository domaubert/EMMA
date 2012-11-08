void dumpmap(int lmap,struct OCT **firstoct,int field,char filename[],REAL zmin, REAL zmax);
void dumpcube(int lmap,struct OCT **firstoct,int field,char filename[],REAL tsim);
void dumppart(struct OCT **firstoct,char filename[],int npart, int levelcoarse, int levelmax,REAL tsim);
void GetParameters(char *fparam, struct RUNPARAMS *param);
struct PART * read_grafic_part(struct PART *part, struct CPUINFO *cpu, REAL *munit, REAL *ainit, REAL *omegam, REAL *omegav, REAL *Hubble, int *npart, REAL omegab);
int read_grafic_hydro(struct CPUINFO *cpu,  REAL *ainit, REAL *omegam, REAL *omegav, REAL *Hubble, REAL *omegab);
void dumpgrid(int levelmax,struct OCT **firstoct, char filename[],REAL tsim);
