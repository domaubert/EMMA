struct PART* findlastpart(struct PART* phead);
int countpart(struct PART* phead);
struct PART* modifpospart(struct PART* phead, REAL len, int dir);
REAL movepart(int levelcoarse,int levelmax,struct OCT** firstoct, REAL dt, struct CPUINFO *cpu);
void  partcellreorg(int levelcoarse,int levelmax,struct OCT **firstoct);
void forcevel(int levelcoarse,int levelmax,struct OCT **firstoct, REAL **vcomp,int stride,REAL dt, struct CPUINFO *cpu, struct PACKET **sendbuffer,struct PACKET **recvbuffer);
REAL egypart(int levelcoarse,int levelmax,struct OCT **firstoct,struct CPUINFO *cpu);
REAL comptstep(int levelcoarse,int levelmax,struct OCT** firstoct, REAL fa, REAL fa2, struct CPUINFO* cpu,REAL);
void accelpart(int level,struct OCT **firstoct, REAL dt, struct CPUINFO *cpu, struct PACKET **sendbuffer,struct PACKET **recvbuffer);
