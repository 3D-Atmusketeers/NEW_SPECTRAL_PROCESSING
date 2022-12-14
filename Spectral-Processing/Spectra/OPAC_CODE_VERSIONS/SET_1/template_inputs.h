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

#define NTEMP 49                           /* Number of temperature points in grid   */
#define NLAMBDA 11215                       /* Number of wavelength points in grid [4616/2598]   */

// This is the Npressure for low res
#define NPRESSURE 28    

#define W0_VAL <<W0_VAL>>
#define G0_VAL <<G0_VAL>>

/* Planet parameters */
#define INPUT_INCLINATION <<inclination>>  /* Planet inclination in radians            */
#define INPUT_PHASE <<phase>>              /* Planet inclination in degrees           */
#define G <<GRAVITY_SI>>                   /* Planet surface gravity                 */

#define R_PLANET <<R_PLANET>>              /* Planet radius at base of atmosphere      */
#define ORB_SEP <<ORB_SEP>>                // This is some distance
#define STELLAR_TEMP <<STELLAR_TEMP>>      // Stellar Blackbody temperature
#define R_STAR <<R_STAR>>                  /* Stellar radius                         */
#define P_ROT  <<P_ROT>>                   /* Rotation period in days (= P_ORB for tidally locked planet)    */
#define HAZE_TYPE <<HAZE_TYPE>>
#define HAZES <<HAZES>>
#define R_VEL 0.0                          /* Radial Velocity                        */
#define MU <<MEAN_MOLECULAR_WEIGHT>>          /* Mean molecular weight                  */
#define FORMAT 2                           /* FORMAT=1 -> small opacity table        */
                                           /* FORMAT=2 -> large opacity table        */

#define ABUND_TRACK_IND 10 /* if test visualize, the index of the mapped species in the EOS file */
#define CHEM_FILE_NCOLS_USED 23 /* number of EOS columns to use (rest discarded). including T and P. */
#define CHEM_FILE_NCOLS 23 /* total number of EOS file columns, including temperature and pressure.*/

/* Opacities for spectra */
#define CHEM_FILE   "DATA/SET_1/eos_solar_GJ1214b.dat"
#define C2H2_FILE   "DATA/SET_1/opacC2H2.dat"
#define CH4_FILE    "DATA/SET_1/opacCH4.dat"
#define CO_FILE     "DATA/SET_1/opacCO.dat"
#define CO2_FILE    "DATA/SET_1/opacCO2.dat"
#define FeH_FILE    "DATA/SET_1/opacFeH.dat"
#define H2O_FILE    "DATA/SET_1/opacH2O.dat"
#define H2S_FILE    "DATA/SET_1/opacH2S.dat"
#define HCN_FILE    "DATA/SET_1/opacHCN.dat"
#define K_FILE      "DATA/SET_1/opacK.dat"
#define Na_FILE     "DATA/SET_1/opacNa.dat"
#define NH3_FILE    "DATA/SET_1/opacNH3.dat"
#define TiO_FILE    "DATA/SET_1/opacTiO.dat"
#define VO_FILE     "DATA/SET_1/opacVO.dat"

#endif /* !__INPUT_H__ */

/* ------- end ---------------------------- input.h  ----------------- */
