#define RANK_DISP (0)
#define FBKP (1)
#define FDUMP (1)
#define GAMMA (5./3.)
#define CFL (0.85)
#define GRAV (0.25)
#define FRACDX  (0.25)

#define SEED (12345)
#define PMIN (1e-12)
#define NCOSMOTAB (128)
#define VBC (30.); // relative DM velocity ala Tseliakovich & Hirata
#define FRAC_VAR (0.1)
#define N_BIN_PDF (128)

#define NGRP_SPACE (1)
#define NGRP_TIME (1)
#define NGRP (NGRP_SPACE * NGRP_TIME)

#define THRESH_Z_XION_MAP (0.5)

// ================= PHYSICAL CONSTANTS ===================
#define LIGHT_SPEED_IN_M_PER_S (299792458.) // m.s-1
#define KBOLTZ (1.3806e-23) // J/K
#define PARSEC (3.085677e16) // m
#define AVOGADRO (6.02214129e23) // mol-1
#define MYR (3.1536e13) // s
#define PROTON_MASS (1.67262158e-27) //kg

#define MHE_OVER_MH (4.002)
#define NEWTON_G (6.67384e-11) // SI
#define HELIUM_MASS (MHE_OVER_MH*PROTON_MASS) //kg
#define SOLAR_MASS (1.989e30) //kg
#define KPLANCK (6.62617e-34)	//J/s
#define ELECTRONVOLT (1.602176565e-19) //J

#ifndef WHYDRO2
  #define OMEGAB (0.0)
#else
  #define OMEGAB (0.046) // 0.049 for PLANCK
#endif

#ifdef HELIUM
#define YHE (0.24) // Helium Mass fraction
#else
#define YHE (0.) // Helium Mass fraction
#endif

#define yHE (YHE/(1.-YHE)/4.002) // Helium number fraction
#define MOLECULAR_MU (1.0)
