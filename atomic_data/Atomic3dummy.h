#define NGRP_ATOMIC 3

#ifdef S_50000
#define SECTION_EFFICACE hnu[0]=20.2761765191*1.6022e-19;alphae[0]=2.49362638578e-22*c;alphai[0]=2.93055200728e-22*c;
#define FACTGRP factgrp[0]=1.0;
#endif

#ifdef S_100000
#define SECTION_EFFICACE hnu[0]=29.609895722*1.6022e-19;alphae[0]=1.09691260769e-22*c;alphai[0]=1.6302904562e-22*c;hnu[1]=29.609895722*1.6022e-19;alphae[1]=1.09691260769e-22*c;alphai[1]=1.6302904562e-22*c;hnu[2]=29.609895722*1.6022e-19;alphae[2]=1.09691260769e-22*c;alphai[2]=1.6302904562e-22*c;
#define FACTGRP factgrp[0]=0.45;factgrp[1]=0.5;factgrp[2]=0.05;
#endif

#ifdef S_1p8
#define SECTION_EFFICACE hnu[0]=29.6259319907*1.6022e-19;alphae[0]=1.47826746468e-22*c;alphai[0]=2.5231304835e-22*c;
#define FACTGRP factgrp[0]=1.0;
#endif

#ifdef S_5p0
#define SECTION_EFFICACE hnu[0]=16.976533961*1.6022e-19;alphae[0]=3.79210609825e-22*c;alphai[0]=4.12837476352e-22*c;
#define FACTGRP factgrp[0]=1.0;
#endif

