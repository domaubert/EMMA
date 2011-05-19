void part2cell_cic(struct PART *curp, struct OCT *curoct, int icell, char full);
void cell2part_cic(struct PART *curp, struct OCT *curoct, int icell, char dir, float dt);
void call_cic(int levelmax,int levelcoarse,struct OCT **firstoct, struct CPUINFO *cpu);
void call_cic2(int levelmax,int levelcoarse,struct OCT **firstoct, struct CPUINFO *cpu);
float cell2part_cic_egy(struct PART *curp, struct OCT *curoct, int icell);
void cell2part_cic_GPU(struct PART *curp, struct OCT *curoct, int icell, char dir, float dt);

  
