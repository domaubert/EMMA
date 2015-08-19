#define RANK_DISP 0
#define FBKP 1
#define GAMMA (5./3.)
#define CFL (0.85)
#define GRAV (0.25)
#define FRACDX  (0.25)
#define PMIN 1e-12
#define NCOSMOTAB (128)
#define VBC (30.); // relative DM velocity ala Tseliakovich & Hirata
#define FRAC_VAR (0.1)
#define SEED (4569) // srand seed for stochastic stars creation

#define NGRP_SPACE (1)
#define NGRP_TIME (1)
#define NGRP (NGRP_SPACE * NGRP_TIME)

#define NVAR_R (5)
#define EMIN (1e-8)
#define YHE (0)



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
#define KPLANCK (6.62617e-34)	//J/s
#define ELECTRONVOLT (1.602176565e-19) //J


// ================= SUPERNOVAE ===========================

// Proportion of ejecta following fig 109 of the Starburst99 model:
// http://www.stsci.edu/science/starburst99/figs/mass_inst_e.html
#define EJECTA_PROP (0.5260172663907063)

// Total feedback energy following fig 115 of the Starburst99 model :
// http://www.stsci.edu/science/starburst99/figs/energy_inst_e.html
#define SN_EGY (9.73565592733335e11) // j/kg



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
