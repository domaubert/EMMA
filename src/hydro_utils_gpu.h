int advancehydroGPU(struct OCT **firstoct, int level, struct CPUINFO *cpu, struct HGRID *stencil, int stride, REAL dxcur, REAL dtnew);
void create_hydstencil_GPU(struct CPUINFO *cpu, int stride);
void create_pinned_stencil(struct HGRID **stencil, int stride);
void destroy_hydstencil_GPU(struct CPUINFO *cpu, int stride);
void destroy_pinned_stencil(struct HGRID **stencil, int stride);
