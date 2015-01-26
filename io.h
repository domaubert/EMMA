void dumpmap(int lmap,struct OCT **firstoct,int field,char filename[],REAL zmin, REAL zmax);
void dumpcube(int lmap,struct OCT **firstoct,int field,char filename[],REAL tsim);
void dumppart(struct OCT **firstoct,char filename[],int levelcoarse, int levelmax,REAL tsim,struct CPUINFO *cpu);
void GetParameters(char *fparam, struct RUNPARAMS *param);
struct PART * read_grafic_part(struct PART *part, struct CPUINFO *cpu, REAL *munit, REAL *ainit, int *npart, struct RUNPARAMS *param, int level);
int read_grafic_hydro(struct CPUINFO *cpu,  REAL *ainit, struct RUNPARAMS *param,int level);
void dumpgrid(int levelmax,struct OCT **firstoct, char filename[],REAL tsim, struct RUNPARAMS *param);
void save_amr(char filename[], struct OCT **firstoct,REAL tsim, REAL, int ,int, struct RUNPARAMS *param, struct CPUINFO *cpu, struct PART *proot, REAL *adt);
struct OCT* restore_amr(char filename[], struct OCT **firstoct, struct OCT **lastoct, REAL *tsim, REAL *,int*,int *,struct RUNPARAMS *param, struct CPUINFO *cpu, struct PART *proot, REAL *adt,struct CELL *root);

void save_part(char filename[],struct OCT **firstoct, int levelcoarse, int levelmax, REAL tsim, struct CPUINFO *cpu, struct PART* proot);
struct PART * restore_part(char filename[], struct OCT **firstoct, REAL *tsim, struct RUNPARAMS *param, struct CPUINFO *cpu, struct PART *proot);
int read_evrard_hydro(struct CPUINFO *cpu,struct OCT **firstoct, struct RUNPARAMS *param);
struct PART * read_edbert_part(struct PART *part, struct CPUINFO *cpu, REAL *munit, REAL *ainit, int *npart, struct RUNPARAMS *param,struct OCT **firstoct);
struct PART * read_zeldovich_part(struct PART *part, struct CPUINFO *cpu, REAL *munit, REAL *ainit, int *npart, struct RUNPARAMS *param, struct OCT **firstoct);
void dumpIO(REAL tsim, struct RUNPARAMS *param,struct CPUINFO *cpu, struct OCT **firstoct, REAL *adt, int pdump);
void read_shocktube(struct CPUINFO *cpu, REAL *ainit, struct RUNPARAMS *param, struct OCT **firstoct);

void dumpHeader(struct RUNPARAMS *param, struct CPUINFO *cpu,char *fparam);

