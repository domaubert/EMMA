#define NCELLS 256
#define NCELLX (NCELLS)
#define NCELLY (NCELLS)
#define NCELLZ (2*NCELLS)
#define NBUFF (NCELLX*NCELLZ*4) // 4 for 4 fields 

#define NGPUX 8
#define NGPUY 8
#define NGPUZ 4

#define BLOCKCOOL 64 //MUST BE <128
#define GRIDCOOLX 16 // can be increased (4 for 256^3)
#define GRIDCOOLY ((NCELLX*NCELLY*NCELLZ)/GRIDCOOLX/BLOCKCOOL) 
