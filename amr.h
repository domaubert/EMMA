struct OCT * refine_cells(int levelcoarse, int levelmax, struct OCT **firstoct, struct OCT ** lastoct, struct OCT * endoct, struct CPUINFO *cpu);
void mark_cells(int levelcoarse,int levelmax,struct OCT **firstoct, int nsmooth, float threshold, struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer);
