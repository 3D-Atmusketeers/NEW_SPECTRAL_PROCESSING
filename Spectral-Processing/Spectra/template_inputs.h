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

#define NTEMP <<num_temperature_points>>                           /* Number of temperature points in grid   */
#define NLAMBDA <<num_wavelength_points>>                       /* Number of wavelength points in grid [4616/2598]   */

// This is the Npressure for low res
#define NPRESSURE <<num_pressure_points>>    /* Number of pressure points in grid   [13/17]   */

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
#define MU 2.36                            /* Mean molecular weight                  */
#define FORMAT 2                           /* FORMAT=1 -> small opacity table        */
                                           /* FORMAT=2 -> large opacity table        */

#define ABUND_TRACK_IND 10 /* if test visualize, the index of the mapped species in the EOS file */
#define CHEM_FILE_NCOLS_USED <<NUM_COLS>> /* number of EOS columns to use (rest discarded). including T and P. */
#define CHEM_FILE_NCOLS <<NUM_COLS>> /* total number of EOS file columns, including temperature and pressure.*/

/* Opacities for spectra */
#define CHEM_FILE   <<CHEM_FILE>>
#define C2H2_FILE   <<C2H2_FILE>>
#define CH4_FILE    <<CH4_FILE>>
#define CO_FILE     <<CO_FILE>>
#define CO2_FILE    <<CO2_FILE>>
#define FeH_FILE    <<FeH_FILE>>
#define H2O_FILE    <<H2O_FILE>>
#define H2S_FILE    <<H2S_FILE>>
#define HCN_FILE    <<HCN_FILE>>
#define K_FILE      <<K_FILE>>
#define Na_FILE     <<Na_FILE>>
#define NH3_FILE    <<NH3_FILE>>
#define TiO_FILE    <<TiO_FILE>>
#define VO_FILE     <<VO_FILE>>

#endif /* !__INPUT_H__ */

/* ------- end ---------------------------- input.h  ----------------- */
