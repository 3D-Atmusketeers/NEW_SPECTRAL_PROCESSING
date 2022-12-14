/*----------------------- input.h ------------------------------------

Defines input values and files for 3-D emission spectra

---------------------------------------------------------------------- */

#ifndef __INPUT_H__
#define __INPUT_H__

/* I/O SETTINGS. */

/* File names */
#define OUTPUT_PREFIX <<output_file>>      /* output name */
#define T_P_3D_FILE <<input_file>>         /* input file */

/* Output settings */
#define N_PHASE 1                          /* Number of phases [96 max; lon grid in increments of 3.75] */
#define DOPPLER <<doppler>>                /* 0:Off; 1:On */
#define CLOUDS <<CLOUDS>>                           /* 0:Off; 1:On */

/* Grid settings */
#define NTAU <<NTAU>>                            /* Number of altitude points in grid      */
#define NLAT  <<NLAT>>                           /* Number of latitude points in 3-D  grid */
#define NLON  <<NLON>>                           /* Number of longitude points in 3-D grid */

#define NTEMP 30                           /* Number of temperature points in grid   */
#define NLAMBDA 2598                       /* Number of wavelength points in grid [4616/2598]   */

// This is the Npressure for low res
#define NPRESSURE 13    /* Number of pressure points in grid   [13/17]   */

#define W0_VAL <<W0_VAL>>
#define G0_VAL <<G0_VAL>>

/* Planet parameters */
#define INPUT_INCLINATION <<inclination>>  /* Planet inclination in radians            */
#define INPUT_PHASE <<phase>>              /* Planet inclination in degrees           */
#define G <<GRAVITY_SI>>                   /* Planet surface gravity                 */

#define R_PLANET <<R_PLANET>>                 /* Planet radius at base of atmosphere      */
#define ORB_SEP <<ORB_SEP>>                  // This is some distance
#define STELLAR_TEMP <<STELLAR_TEMP>>                // Stellar Blackbody temperature
#define R_STAR <<R_STAR>>                    /* Stellar radius                         */
#define P_ROT  <<P_ROT>>                        /* Rotation period in days (= P_ORB for tidally locked planet)    */
#define HAZE_TYPE <<HAZE_TYPE>>
#define HAZES <<HAZES>>

#define R_VEL 0.0                          /* Radial Velocity                        */
#define MU <<MEAN_MOLECULAR_WEIGHT>>          /* Mean molecular weight                  */
#define FORMAT 2                           /* FORMAT=1 -> small opacity table        */
                                           /* FORMAT=2 -> large opacity table        */

/* Aerosol properties (calculated by the Mischenko Mie code) */

#define PI0_KCl 0.74
#define G0_KCl 0.15
#define QE_KCl 0.12
#define PI0_ZnS 0.74
#define G0_ZnS 0.15
#define QE_ZnS 0.12
#define PI0_Na2S 0.74
#define G0_Na2S 0.15
#define QE_Na2S 0.12
#define PI0_MnS 0.74
#define G0_MnS 0.15
#define QE_MnS 0.12
#define PI0_Cr 0.74
#define G0_Cr 0.15
#define QE_Cr 0.12
#define PI0_SiO2 0.74
#define G0_SiO2 0.15
#define QE_SiO2 0.12
#define PI0_Mg2SiO4 0.74
#define G0_Mg2SiO4 0.15
#define QE_Mg2SiO4 0.12
#define PI0_VO 0.74
#define G0_VO 0.15
#define QE_VO 0.12
#define PI0_Ni 0.74
#define G0_Ni 0.15
#define QE_Ni 0.12
#define PI0_Fe 0.74
#define G0_Fe 0.15
#define QE_Fe 0.12
#define PI0_CaSiO4 0.74
#define G0_CaSiO4 0.15
#define QE_CaSiO4 0.12
#define PI0_CaTiO3 0.74
#define G0_CaTiO3 0.15
#define QE_CaTiO3 0.15
#define PI0_Al2O3 0.74
#define G0_Al2O3 0.15
#define QE_Al2O3 0.12

/* Opacities for spectra */
#define CHEM_FILE   "DATA/SET_2/eos_solar_doppler.dat"
#define CH4_FILE    "DATA/SET_2/opacCH4.dat"
#define CO2_FILE    "DATA/SET_2/opacCO2.dat"
#define CO_FILE     "DATA/SET_2/opacCO.dat"
#define H2O_FILE    "DATA/SET_2/opacH2O.dat"
#define NH3_FILE    "DATA/SET_2/opacNH3.dat"
#define O2_FILE     "DATA/SET_2/opacO2.dat"
#define O3_FILE     "DATA/SET_2/opacO3.dat"

#endif /* !__INPUT_H__ */

/* ------- end ---------------------------- input.h  ----------------- */