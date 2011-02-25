
int segment_cell(struct OCT *curoct, int icell, struct CPUINFO *cpu, int levelcoarse);
void assigncpu2coarseoct(struct OCT *curoct, struct CPUINFO *cpu, int levelcoarse);
int segment_part(float xc,float yc,float zc, struct CPUINFO *cpu, int levelcoarse);
int hfun(unsigned key);
void load_balance(int levelcoarse,struct CPUINFO *cpu);

