struct PART* findlastpart(struct PART* phead);
int countpart(struct PART* phead);
struct PART* modifpospart(struct PART* phead, float len, int dir);
float movepart(int levelcoarse,int levelmax,struct OCT** firstoct, float dt, struct CPUINFO *cpu);
void  partcellreorg(int levelcoarse,int levelmax,struct OCT **firstoct);
void forcevel(int levelcoarse,int levelmax,struct OCT **firstoct, float **vcomp,int stride,float dt, struct CPUINFO *cpu, struct PACKET **sendbuffer,struct PACKET **recvbuffer);
float egypart(int levelcoarse,int levelmax,struct OCT **firstoct,struct CPUINFO *cpu);
float comptstep(int levelcoarse,int levelmax,struct OCT** firstoct, float dt, struct CPUINFO* cpu);
