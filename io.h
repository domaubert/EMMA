void dumpmap(int lmap,struct OCT **firstoct,int field,char filename[],float zmin, float zmax);
void dumpcube(int lmap,struct OCT **firstoct,int field,char filename[],float tsim);
void dumppart(struct OCT **firstoct,char filename[],int npart, int levelcoarse, int levelmax,float tsim);
void GetParameters(char *fparam, struct RUNPARAMS *param);
