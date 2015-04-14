#define RANK_DISP 0
#define FBKP 1
#define GAMMA (5./3.)
#define CFL (0.85)
#define GRAV (0.25)
#define FRACDX  (0.25)
#define SEED (12345)
#define PMIN 1e-12
#define NCOSMOTAB (262144)
#define VBC (30.); // relative DM velocity ala Tseliakovich & Hirata
#define FRAC_VAR (0.1)
#define SEED (4569) // srand seed for stochastic stars creation
#define NGRP (1)
#define NVAR_R (5)
#define EMIN (1e-8)




// ================= PHYSICAL CONSTANTS ===================
#define LIGHT_SPEED_IN_M_PER_S (299792458.)
#define KBOLTZ (1.3806e-23) // J/K
#define PARSEC (3.085677e16) // in m
#define AVOGADRO (6.02214129e23) // mol-1
#define MYR (3.1536e13) // s
#define PROTON_MASS (1.67262158e-27) //kg
#define NEWTON_G (6.67384e-11) // SI
#define HELIUM_MASSFRACTION (0.24)
#define MOLECULAR_MU (1.0)
#define SOLAR_MASS (1.989e30) //kg


// ================= SUPERNOVAE ===========================

//#define SN_EGY (3.7e11) 		// 3.7e15 erg.g-1 -> 3.7e11 J.kg-1 ->  Kay 2002   // 4e48 erg.Mo-1 springel hernquist 2003 -> OK
//#define N_SNII (1.)

//Vecchia & Schaye 2012
//efine SN_EGY (8.73e11/1.736e-2)
//chabrier IMF
//efine N_SNII (1.736e-2)    //[6->100MO]
//#define N_SNII (1.180e-2)    //[8->100MO]

//Salpeter IMF
//#define N_SNII (1.107e-2)   //[6->100MO]
//#define N_SNII (0.742e-2)    //[8->100MO]