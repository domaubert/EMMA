struct OCT* cell2oct(struct CELL* cell);
void getcellnei(int cindex, int *neip, int *cell);
void cic_child(struct OCT* oct,struct OCT* octorg, int icellorg);
void setOctList(struct OCT *firstoct, struct CPUINFO *cpu, struct RUNPARAMS *param, int level);
void oct2loct(struct OCT *oct, struct LOCT *loct);
