
void remove_valvec_GPU(float *vec, int nval, int stride, float avg, int level, int *vecl);
float square_vec_GPU(float *vec, int nval, int stride, int level, int *vecl, float *vec2, float *vecsum);
float laplacian_vec2_GPU(float *vecden,float *vecpot,float *vecpotnew,int *vecnei,int *vecl, int *vecicoarse, int level, int nread,int stride,float dx,float omegam,float tsim, float *vres, float *vecsum);