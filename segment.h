
int segment_cell(struct OCT *curoct, int icell, struct CPUINFO *cpu, int levelcoarse);
void assigncpu2coarseoct(struct OCT *curoct, struct CPUINFO *cpu, int levelcoarse);
int segment_part(REAL xc,REAL yc,REAL zc, struct CPUINFO *cpu, int levelcoarse);
unsigned long hfun(unsigned long key, unsigned long maxval);
void load_balance(int levelcoarse,struct CPUINFO *cpu);
unsigned long oct2key(struct OCT *curoct,int level);
unsigned long pos2key(REAL xc, REAL yc, REAL zc, int level);
