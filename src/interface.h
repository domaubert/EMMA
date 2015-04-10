void CPU2GPU(REAL *gpupt, REAL *cpupt, int noctet);
void GPU2CPU(REAL *gpupt, REAL *cpupt, int noctet);
void GPU2GPU(REAL *gpupt, REAL *cpupt, int noctet);
void CPU2GPU_INT(int *gpupt, int *cpupt, int noctet);
void GPU2CPU_INT(int *gpupt, int *cpupt, int noctet);
void CPU2GPU_UINT(unsigned int *gpupt, unsigned int *cpupt, int noctet);
void GPU2CPU_UINT(unsigned int *gpupt, unsigned int *cpupt, int noctet);
REAL * GPUallocREAL(int nmem);
unsigned long int GPUallocScanPlan(int stride);
