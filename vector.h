
struct OCT* gathercomp(struct OCT *octstart, float *vec, char nei, char var, int stride, struct CPUINFO *cpu, int *nread);
struct OCT* scattercomp(struct OCT *octstart, float *vec, char nei, char var, int stride, struct CPUINFO *cpu);
float laplacian(float **vcomp, int stride, float dx, int locres);
float laplaciancosmo(float **vcomp, int stride, float dx, int locres, float omegam, float a);
void grad(float **vcomp, int stride, float dx, int dir);
void remove_avg(float *vcomp, int stride, float avg);
float square(float *vcomp, int stride);
float square_res(float **vcomp, int stride, float dx);

float laplacian_vec(float *vecden,float *vecpot,float *vecpotnew,int *vecnei, int *vecl, int level, int nread,int stride,float dx,float omegam,float tsim);
float laplacian_vec2(float *vecden,float *vecpot,float *vecpotnew,int *vecnei, int *vecl, int*, int *,int level, int curcpu, int nread,int stride,float dx,float factdens);
float square_vec(float *vec, int nval, int stride, int, int, int*, int *);
void remove_valvec(float *vec, int nval, int stride, float avg, int, int *);
struct OCT *gathervec(struct OCT *octstart, float *vec, char var, int *vecl, int stride, struct CPUINFO *cpu, int *nread);
struct OCT *gathervec2(struct OCT *octstart, float *vec, char var, int *vecl, int *,int *veccpu, int stride, struct CPUINFO *cpu, int *nread);
struct OCT *gathervec2_light(struct OCT *octstart, float *vec, char var, int stride, struct CPUINFO *cpu, int *nread, int level);

struct OCT *gathervecnei(struct OCT *octstart, int *vecnei, float *vec, char var, int *vecl, int stride, struct CPUINFO *cpu, int *nread);
struct OCT *gathervecnei2(struct OCT *octstart, int *vecnei, int stride, struct CPUINFO *cpu, int *nread);
int countvecocts(struct OCT *octstart, int stride, struct CPUINFO *cpu, int *nread);

struct OCT *checknei(struct OCT *octstart, int *vecnei, int stride);
struct OCT *scattervec(struct OCT *octstart, float *vec, char var, int stride, struct CPUINFO *cpu, int nread);
struct OCT *scattervec_light(struct OCT *octstart, float *vec, char var, int stride, struct CPUINFO *cpu, int nread,int level);
int residual_vec2(float *vecden,float *vecpot,float *vecres,int *vecnei,int *vecl, int *vecicoarse, int *veccpu, int level, int curcpu, int nread,int stride,float dx,float factdens);
struct OCT *gathervechydro(struct OCT *octstart, struct MULTIVECT *data, int stride, struct CPUINFO *cpu, int *nread);
struct OCT *scattervechydro(struct OCT *octstart, struct MULTIVECT *data, int stride, struct CPUINFO *cpu);
