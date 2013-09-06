
#include <cudpp.h>

#ifdef WCUDA_ERR
#  define CUDA_CHECK_ERROR(errorMessage) {                                    \
    cudaError_t err = cudaGetLastError();                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
                errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    }
#else
#define CUDA_CHECK_ERROR(errorMessage)
#endif

// ====================== structure for CUDPP =======
struct CUPARAM{
  CUDPPHandle theCudpp;
  CUDPPConfiguration config;
  CUDPPHandle scanplan;
};
