
__global__ void cusetboundaryref_xp(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx);
__global__ void cusetboundaryref_yp(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx);
__global__ void cusetboundaryref_zp(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx);
__global__ void cusetboundaryref_xm(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx);
__global__ void cusetboundaryref_ym(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx);
__global__ void cusetboundaryref_zm(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx);

__global__ void cusetboundarytrans_xp(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx);
__global__ void cusetboundarytrans_yp(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx);
__global__ void cusetboundarytrans_zp(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx);
__global__ void cusetboundarytrans_xm(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx);
__global__ void cusetboundarytrans_ym(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx);
__global__ void cusetboundarytrans_zm(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx);

__global__ void cusetboundaryper_xp(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx);
__global__ void cusetboundaryper_yp(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx);
__global__ void cusetboundaryper_zp(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx);
__global__ void cusetboundaryper_xm(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx);
__global__ void cusetboundaryper_ym(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx);
__global__ void cusetboundaryper_zm(float *cuegy, float *cuxion, float *cudensity, float *cutemperature, float *cuflx);
