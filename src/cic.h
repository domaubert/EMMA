void part2cell_cic(struct PART *curp, struct OCT *curoct, int icell, char full);
void cell2part_cic(struct PART *curp, struct OCT *curoct, int icell, REAL dt);
void call_cic(int levelmax,int levelcoarse,struct OCT **firstoct, struct CPUINFO *cpu);
void call_cic2(int levelmax,int levelcoarse,struct OCT **firstoct, struct CPUINFO *cpu);
REAL cell2part_cic_egy(struct PART *curp, struct OCT *curoct, int icell);
void cell2part_cic_GPU(struct PART *curp, struct OCT *curoct, int icell, char dir, REAL dt);
void L_clean_dens(int level,struct RUNPARAMS *param, struct OCT **firstoct, struct CPUINFO *cpu);
void L_cic(int level,struct OCT **firstoct, struct RUNPARAMS *param, struct CPUINFO *cpu);

