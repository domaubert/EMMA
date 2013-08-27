void dumpmap(int lmap,struct OCT **firstoct,int field,char filename[],REAL zmin, REAL zmax);
void dumpcube(int lmap,struct OCT **firstoct,int field,char filename[],REAL tsim);
void dumppart(struct OCT **firstoct,char filename[],int levelcoarse, int levelmax,REAL tsim,struct CPUINFO *cpu);
void GetParameters(char *fparam, struct RUNPARAMS *param);
struct PART * read_grafic_part(struct PART *part, struct CPUINFO *cpu, REAL *munit, REAL *ainit, int *npart, struct RUNPARAMS *param);
int read_grafic_hydro(struct CPUINFO *cpu,  REAL *ainit, struct RUNPARAMS *param);
void dumpgrid(int levelmax,struct OCT **firstoct, char filename[],REAL tsim, struct RUNPARAMS *param);
void save_amr(char filename[], struct OCT **firstoct,REAL tsim, REAL, int ,int, struct RUNPARAMS *param, struct CPUINFO *cpu, struct PART *proot, REAL *adt);
struct OCT* restore_amr(char filename[], struct OCT **firstoct, struct OCT **lastoct, REAL *tsim, REAL *,int*,int *,struct RUNPARAMS *param, struct CPUINFO *cpu, struct PART *proot, REAL *adt);

void save_part(char filename[],struct OCT **firstoct, int levelcoarse, int levelmax, REAL tsim, struct CPUINFO *cpu, struct PART* proot);
struct PART * restore_part(char filename[], struct OCT **firstoct, REAL *tsim, struct RUNPARAMS *param, struct CPUINFO *cpu, struct PART *proot);
