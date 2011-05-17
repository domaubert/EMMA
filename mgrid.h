
float poisson_jacob(int level,int levelcoarse,int levelmax, struct OCT **firstoct,struct MULTIVECT* vectors,int stride,struct CPUINFO *cpu,float omegam,float tsim, struct PACKET **sendbuffer,struct PACKET **recvbuffer, int niter, float acc);
float  poisson_mgrid(int level,int levelcoarse,int levelmax,int levelmin, struct OCT **firstoct,struct MULTIVECT *vectors,int stride, struct CPUINFO *cpu, float omegam, float tsim, struct PACKET** sendbuffer, struct PACKET **recvbuffer, int niter, float acc);
