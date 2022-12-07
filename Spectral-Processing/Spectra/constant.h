/*----------------------- constant.h --------------------------------

Author: Eliza Miller-Ricci (emillerricci@cfa.harvard.edu)
Last modified: February 21, 2008

--------------------------------------------------------------------- */

#ifndef __CONSTANT_H__
#define __CONSTANT_H__

/*---- Physical constants ------------------------------------------- */

#define CLIGHT      2.99792458E+08 /* Speed of light [m/s]            */
#define HPLANCK     6.6260755E-34  /* Planck's constant [Js]          */
#define KBOLTZMANN  1.380658E-23   /* Boltzmann's constant[J/K]       */
#define MELECTRON   9.1093897E-31  /* Electron mass [kg]              */
#define QELECTRON   1.60217733E-19 /* Electron charge [C]             */
#define AMU         1.6605402E-27  /* Atomic mass unit [kg]           */
#define MHATOM      1.673725E-27   /* H atom mass [kg]                */
#define NAVOGADRO   6.0221420E+23  /* Avogadro's number [1/mole]      */
#define R_GAS       8.314472E+00   /* Ideal gas const. [J/mol/K]      */
#define SIGMA       5.6704E-08     /* Stef.-Boltz. const. [W/m^2/K^4] */
#define NL          2.6867774E+25  /* Loschmidt constant in [m^-3]    */

/*---- Unit conversions --------------------------------------------- */

#define NM_TO_M       1.0E-09
#define CM_TO_M       1.0E-02
#define KM_TO_M       1.0E+03
#define MICRON_TO_NM  1.0E+03 
#define ERG_TO_JOULE  1.0E-07
#define CAL_TO_JOULE  4.1840
#define G_TO_KG       1.0E-03
#define ATM_TO_PA     1.01325E+05

/*---- Mathematical constants --------------------------------------- */

#ifndef PI
#define PI          3.14159265358979
#endif
#define SQRTPI      1.77245385090551

#endif /* !__CONSTANT_H__ */

/*---- end ------------------------ constant.h ---------------------- */
