
struct OCT* gathercomp(struct OCT *octstart, REAL *vec, char nei, char var, int stride, struct CPUINFO *cpu, int *nread);
struct OCT* scattercomp(struct OCT *octstart, REAL *vec, char nei, char var, int stride, struct CPUINFO *cpu);
REAL laplacian(REAL **vcomp, int stride, REAL dx, int locres);
REAL laplaciancosmo(REAL **vcomp, int stride, REAL dx, int locres, REAL omegam, REAL a);
void grad(REAL **vcomp, int stride, REAL dx, int dir);
void remove_avg(REAL *vcomp, int stride, REAL avg);
REAL square(REAL *vcomp, int stride);
REAL square_res(REAL **vcomp, int stride, REAL dx);

REAL laplacian_vec(REAL *vecden,REAL *vecpot,REAL *vecpotnew,int *vecnei, int *vecl, int level, int nread,int stride,REAL dx,REAL omegam,REAL tsim);
REAL laplacian_vec2(REAL *vecden,REAL *vecpot,REAL *vecpotnew,int *vecnei, int *vecl, int*, int *,int level, int curcpu, int nread,int stride,REAL dx,REAL factdens);
REAL square_vec(REAL *vec, int nval, int stride, int, int, int*, int *);
void remove_valvec(REAL *vec, int nval, int stride, REAL avg, int, int *);
struct OCT *gathervec(struct OCT *octstart, REAL *vec, char var, int *vecl, int stride, struct CPUINFO *cpu, int *nread);
struct OCT *gathervec2(struct OCT *octstart, REAL *vec, char var, int *vecl, int *,int *veccpu, int stride, struct CPUINFO *cpu, int *nread);
struct OCT *gathervec2_light(struct OCT *octstart, REAL *vec, char var, int stride, struct CPUINFO *cpu, int *nread, int level);

struct OCT *gathervecnei(struct OCT *octstart, int *vecnei, REAL *vec, char var, int *vecl, int stride, struct CPUINFO *cpu, int *nread);
struct OCT *gathervecnei2(struct OCT *octstart, int *vecnei, int stride, struct CPUINFO *cpu, int *nread);
int countvecocts(struct OCT *octstart, int stride, struct CPUINFO *cpu, int *nread);

struct OCT *checknei(struct OCT *octstart, int *vecnei, int stride);
struct OCT *scattervec(struct OCT *octstart, REAL *vec, char var, int stride, struct CPUINFO *cpu, int nread);
struct OCT *scattervec_light(struct OCT *octstart, REAL *vec, char var, int stride, struct CPUINFO *cpu, int nread,int level);
int residual_vec2(REAL *vecden,REAL *vecpot,REAL *vecres,int *vecnei,int *vecl, int *vecicoarse, int *veccpu, int level, int curcpu, int nread,int stride,REAL dx,REAL factdens);
struct OCT *gathervechydro(struct OCT *octstart, struct MULTIVECT *data, int stride, struct CPUINFO *cpu, int *nread);
struct OCT *scattervechydro(struct OCT *octstart, struct MULTIVECT *data, int stride, struct CPUINFO *cpu);
struct OCT *gatherstencil(struct OCT *octstart, struct HGRID *stencil, int stride, struct CPUINFO *cpu, int *nread);
struct OCT *scatterstencil(struct OCT *octstart, struct HGRID *stencil, int stride, struct CPUINFO *cpu, REAL dxcur, REAL dtnew);

