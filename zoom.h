void zoom_level(int level, struct CPUINFO *cpu, struct RUNPARAMS *param, struct OCT **firstoct,  struct OCT ** lastoct);
int queryzoom(struct OCT *curoct, int icell, REAL dxcur, REAL Rin);
int pos2levelzoom(REAL xc, REAL yc, REAL zc, struct RUNPARAMS *param);
