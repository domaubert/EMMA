#define RANK_DISP 0
#define FBKP 1
#define GAMMA (5./3.)
#define CFL (0.85)
#define GRAV (0.25)
#define FRACDX  (0.25)

#define SEED (12345)
#define PMIN 1e-9
#define NCOSMOTAB (128)
#define VBC (30.); // relative DM velocity ala Tseliakovich & Hirata
#define FRAC_VAR (0.1)


// ================= PHYSICAL CONSTANTS ===================
#define LIGHT_SPEED_IN_M_PER_S (299792458.)
#define KBOLTZ (1.3806e-23) // J/K
#define PARSEC (3.085677e16) // in m
#define AVOGADRO (6.02214129e23) // mol-1
#define MYR (3.1536e13) // s
#define PROTON_MASS (1.67262158e-27) //kg

#define MHE_OVER_MH (4.002)
#define HELIUM_MASS (MHE_OVER_MH*PROTON_MASS) //kg
#define NEWTON_G (6.67384e-11) // SI

#ifndef WHYDRO2
  #define OMEGAB (0.0)
#else
  #define OMEGAB (0.049) // 0.049 for PLANCK
#endif

#ifdef HELIUM
#define YHE (0.24) // Helium Mass fraction
#else
#define YHE (0.) // Helium Mass fraction
#endif

#define yHE (YHE/(1.-YHE)/4.002) // Helium number fraction
#define MOLECULAR_MU (1.0)
#define SOLAR_MASS (1.989e30) //kg
#define KPLANCK (6.62617e-34)	//J/s
#define ELECTRONVOLT (1.602176565e-19) //J
#define HELIUM_MASSFRACTION (0.24) // should be removed ?


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
