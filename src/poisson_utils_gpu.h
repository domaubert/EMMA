REAL PoissonJacobiGPU(int level,struct RUNPARAMS *param, struct OCT ** firstoct,  struct CPUINFO *cpu, struct STENGRAV *stencil, int stride, REAL tsim);
void create_gravstencil_GPU(struct CPUINFO *cpu, int stride);
//void create_pinned_gravstencil(struct STENGRAV *gstencil, int stride);
void destroy_pinned_gravstencil(struct STENGRAV *gstencil, int stride);
void destroy_gravstencil_GPU(struct CPUINFO *cpu, int stride);
