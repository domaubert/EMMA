
struct OCT* gathercomp(struct OCT *octstart, float *vec, char nei, char var, int stride, struct CPUINFO *cpu, int *nread);
struct OCT* scattercomp(struct OCT *octstart, float *vec, char nei, char var, int stride, struct CPUINFO *cpu);
void laplacian(float **vcomp, int stride, float dx);
void grad(float **vcomp, int stride, float dx, int dir);
void remove_avg(float *vcomp, int stride, float avg);
float square(float *vcomp, int stride);
float square_res(float **vcomp, int stride, float dx);
