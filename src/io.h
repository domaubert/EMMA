void dumpIO(REAL tsim, struct RUNPARAMS *param,struct CPUINFO *cpu, struct OCT **firstoct, REAL *adt, int pdump);

float assign_grid_field(int field,struct CELL *cell);

void create_HDF5_file(struct CPUINFO *cpu);


void set_offset(struct RUNPARAMS *param, struct CPUINFO *cpu);
void dump_domain(struct RUNPARAMS *param, struct CPUINFO *cpu);
