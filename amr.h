struct OCT * refine_cells(int levelcoarse, int levelmax, struct OCT **firstoct, struct OCT ** lastoct, struct OCT * endoct, struct CPUINFO *cpu,struct OCT * );
void mark_cells(int levelcoarse,int levelmax,struct OCT **firstoct, int nsmooth, REAL threshold, struct CPUINFO *cpu, struct PACKET **sendbuffer, struct PACKET **recvbuffer);
void clean_marks(int levelmax,struct OCT **firstoct);
