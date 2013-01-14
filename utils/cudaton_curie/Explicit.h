
__global__ void cuComputeIon(float *cuegy_new, float *cuflx_new, float *cuxion, float *cudensity, float *cutemperature, float dt, float c, float egy_min, float unit_number, float aexp);
__global__ void cuComputeTemp(float *cuxion, float *cudensity, float *cutemperature, float *cuegy_new, float fudgecool, float c, float dt,float unit_number, int ncvgcool,float aexp,float hubblet, float *,float, float, float,float,float *);

//**********************************************************
#ifdef SDISCRETE
__global__ void cuAddSource(float *cuegy_new,float *cuflx_new, float *cusrc0, int *cusrc0pos,float dt, float dx, int nsource,float aexp,float c);
__global__ void cuSpotSource(float *cuxion, int *cusrc0pos, int mainsource);
#endif

__global__ void cuComputeELF(float *cuegy, float *cuflx, float *cuegy_new, float c, float dx, float dt, int iter, float aexp,float egy_min);
__global__ void cuComputeF_TOTAL_LF(float *cuflx, float *cuflx_new, float c, float dx, float dt, int iter, float *cuegy, float aexp);
