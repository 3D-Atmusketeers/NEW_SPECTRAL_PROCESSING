/*----------------------- input.h ------------------------------------

Defines input values and files for 3-D emission spectra

---------------------------------------------------------------------- */

#ifndef __INPUT_H__
#define __INPUT_H__

/* I/O SETTINGS. */

/* File names */
#define OUTPUT_PREFIX "OUT/Spec_0_HD209_PICKET_NUC_CLOUDS_phase_0.0_inc_0.00.00.0"      /* output name */
#define T_P_3D_FILE "DATA/init_HD209_PICKET_NUC_CLOUDS_phase_0.0_inc_0.0.txt"         /* input file */

/* Output settings */
#define N_PHASE 1                          /* Number of phases [96 max; lon grid in increments of 3.75] */
#define DOPPLER 0                /* 0:Off; 1:On */
#define CLOUDS 1                           /* 0:Off; 1:On */

/* Grid settings */
#define NTAU 250                            /* Number of altitude points in grid      */
#define NLAT  48                           /* Number of latitude points in 3-D  grid */
#define NLON  96                           /* Number of longitude points in 3-D grid */

#define NTEMP 30                           /* Number of temperature points in grid   */
#define NLAMBDA 2598                       /* Number of wavelength points in grid [4616/2598]   */

// This is the Npressure for low res
#define NPRESSURE 17    /* Number of pressure points in grid   [13/17]   */

#define W0_VAL 0.0
#define G0_VAL 0.0

/* Planet parameters */
#define INPUT_INCLINATION 0.0  /* Planet inclination in radians            */
#define INPUT_PHASE 0.0              /* Planet inclination in degrees           */
#define G 8.979999542236328                   /* Planet surface gravity                 */

#define R_PLANET 98531992.0                 /* Planet radius at base of atmosphere      */
#define ORB_SEP 7293758324.342423                  // This is some distance
#define STELLAR_TEMP 6091.0                // Stellar Blackbody temperature
#define R_STAR 827883000.0                    /* Stellar radius                         */
#define P_ROT  3.525001717240822                        /* Rotation period in days (= P_ORB for tidally locked planet)    */
#define HAZE_TYPE "None"
#define HAZES 0

#define R_VEL 0.0                          /* Radial Velocity                        */
#define MU 2.3601          /* Mean molecular weight                  */
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
#define CHEM_FILE   "DATA/SET_3/eos_solar_doppler.dat"
#define CH4_FILE    "DATA/SET_3/opacCH4_hires.dat"
#define CO2_FILE    "DATA/SET_3/opacCO2_hires.dat"
#define CO_FILE     "DATA/SET_3/opacCO_hires.dat"
#define H2O_FILE    "DATA/SET_3/opacH2O_hires.dat"
#define NH3_FILE    "DATA/SET_3/opacNH3_hires.dat"
#define O2_FILE     "DATA/SET_3/opacO2_hires.dat"
#define O3_FILE     "DATA/SET_3/opacO3_hires.dat"

#endif /* !__INPUT_H__ */

/* ------- end ---------------------------- input.h  ----------------- */