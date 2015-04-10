int advanceradGPU (struct OCT **firstoct, int level, struct CPUINFO *cpu, struct RGRID *stencil, int stride, REAL dxcur, REAL dtnew,REAL aexp, struct RUNPARAMS *param, int chemonly);
void create_radstencil_GPU(struct CPUINFO *cpu, int stride);
void create_pinned_stencil_rad(struct RGRID **stencil, int stride);
void destroy_radstencil_GPU(struct CPUINFO *cpu, int stride);
void destroy_pinned_stencil_rad(struct RGRID **stencil, int stride);
void create_param_GPU(struct RUNPARAMS *param, struct CPUINFO *cpu);
