
struct OCT* gathercomp(struct OCT *octstart, float *vec, char nei, char var, int stride, struct CPUINFO *cpu, int *nread);
struct OCT* scattercomp(struct OCT *octstart, float *vec, char nei, char var, int stride, struct CPUINFO *cpu);
float laplacian(float **vcomp, int stride, float dx, int locres);
float laplaciancosmo(float **vcomp, int stride, float dx, int locres, float omegam, float a);
void grad(float **vcomp, int stride, float dx, int dir);
void remove_avg(float *vcomp, int stride, float avg);
float square(float *vcomp, int stride);
float square_res(float **vcomp, int stride, float dx);

float laplacian_vec(float *vecden,float *vecpot,float *vecpotnew,int *vecnei,int nread,int stride,float dx,float omegam,float tsim);
float square_vec(float *vec, int nval, int stride);
void remove_valvec(float *vec, int nval, int stride, float avg);
struct OCT *gathervec(struct OCT *octstart, float *vec, char var, int stride, struct CPUINFO *cpu, int *nread);
struct OCT *gathervecnei(struct OCT *octstart, int *vecnei, int stride, struct CPUINFO *cpu, int *nread);
struct OCT *checknei(struct OCT *octstart, int *vecnei, int stride);
struct OCT *scattervec(struct OCT *octstart, float *vec, char var, int stride, struct CPUINFO *cpu, int nread);
