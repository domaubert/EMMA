int hydrosolve(struct MULTIVECT *data, int level, int curcpu, int nread,int stride,float dx, float dt);
float comptstep_hydro(int levelcoarse,int levelmax,struct OCT** firstoct, float fa, float fa2, struct CPUINFO* cpu, float tmax);
