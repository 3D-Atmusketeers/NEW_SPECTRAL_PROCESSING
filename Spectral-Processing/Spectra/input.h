/*----------------------- input.h ------------------------------------

Defines input values and files for 3-D emission spectra

---------------------------------------------------------------------- */

#ifndef __INPUT_H__
#define __INPUT_H__

/* I/O SETTINGS. */

/* File names */
#define OUTPUT_PREFIX "OUT/Spec_0_GJ1214b-none-0clouds-1met_phase_0.0_inc_0.0"      /* output name */
#define T_P_3D_FILE "DATA/init_GJ1214b-none-0clouds-1met_phase_0.0_inc_0.0.txt"         /* input file */

/* Output settings */
#define N_PHASE 1                          /* Number of phases [96 max; lon grid in increments of 3.75] */
#define DOPPLER 0                /* 0:Off; 1:On */
#define CLOUDS 0                           /* 0:Off; 1:On */

/* Grid settings */
#define NTAU 250                            /* Number of altitude points in grid      */
#define NLAT  48                           /* Number of latitude points in 3-D  grid */
#define NLON  96                           /* Number of longitude points in 3-D grid */

#define NTEMP 49                              /* Number of temperature points in grid   */
#define NLAMBDA 36891                   /* Number of wavelength points in grid [4616/2598]   */
#define LAMBDA_START 13864       /*The start of the lambda index you want to calculate*/
#define LAMBDA_END 34015           /*The end of the lambda index you want to calculate*/

// This is the Npressure for low res
#define NPRESSURE 28    

#define W0_VAL 0.0
#define G0_VAL 0.0

/* Planet parameters */
#define INPUT_INCLINATION 0.0  /* Planet inclination in radians            */
#define INPUT_PHASE 0.0              /* Planet inclination in degrees           */
#define G 10.65                   /* Planet surface gravity                 */

#define R_PLANET 17469282.0              /* Planet radius at base of atmosphere      */
#define ORB_SEP 2137784000.0                // This is some distance
#define STELLAR_TEMP 3250.0      // Stellar Blackbody temperature
#define R_STAR 144009900.0                  /* Stellar radius                         */
#define P_ROT  1.57577577825418                   /* Rotation period in days (= P_ORB for tidally locked planet)    */
#define HAZE_TYPE "None"
#define HAZES 0
#define R_VEL 0.0                          /* Radial Velocity                        */
#define MU 2.3267          /* Mean molecular weight                  */
#define FORMAT 2                           /* FORMAT=1 -> small opacity table        */
                                           /* FORMAT=2 -> large opacity table        */

#define ABUND_TRACK_IND 2 /* if test visualize, the index of the mapped species in the EOS file */
#define CHEM_FILE_NCOLS 47

/* Opacities for spectra */
#define CHEM_FILE   "DATA/chemistry_grid/fastchem_grid_allspecies_ions_lotemp_Z_0.0_C_to_O_1.0.dat"

#define CIA_FILE "DATA/Low-Res/opacCIA.dat"
#define Ca_plus_FILE   "DATA/Low-Res/opacCa_plus.dat"
#define Mg_FILE   "DATA/Low-Res/opacMg.dat"
#define SiH4_FILE   "DATA/Low-Res/opacSiH4.dat"
#define Fe_plus_FILE   "DATA/Low-Res/opacFe_plus.dat"
#define NaH_FILE   "DATA/Low-Res/opacNaH.dat"
#define O2_FILE   "DATA/Low-Res/opacO2.dat"
#define PH3_FILE   "DATA/Low-Res/opacPH3.dat"
#define TiO_FILE   "DATA/Low-Res/opacTiO.dat"
#define CaO_FILE   "DATA/Low-Res/opacCaO.dat"
#define C_FILE   "DATA/Low-Res/opacC.dat"
#define OH_FILE   "DATA/Low-Res/opacOH.dat"
#define CrH_FILE   "DATA/Low-Res/opacCrH.dat"
#define SiH_FILE   "DATA/Low-Res/opacSiH.dat"
#define SO2_FILE   "DATA/Low-Res/opacSO2.dat"
#define CO_FILE   "DATA/Low-Res/opacCO.dat"
#define HF_FILE   "DATA/Low-Res/opacHF.dat"
#define H2S_FILE   "DATA/Low-Res/opacH2S.dat"
#define NH3_FILE   "DATA/Low-Res/opacNH3.dat"
#define NO_FILE   "DATA/Low-Res/opacNO.dat"
#define CH4_FILE   "DATA/Low-Res/opacCH4.dat"
#define SiO_FILE   "DATA/Low-Res/opacSiO.dat"
#define HCN_FILE   "DATA/Low-Res/opacHCN.dat"
#define MgH_FILE   "DATA/Low-Res/opacMgH.dat"
#define SH_FILE   "DATA/Low-Res/opacSH.dat"
#define CO2_FILE   "DATA/Low-Res/opacCO2.dat"
#define VO_FILE   "DATA/Low-Res/opacVO.dat"
#define H2O_FILE   "DATA/Low-Res/opacH2O.dat"
#define AlH_FILE   "DATA/Low-Res/opacAlH.dat"
#define CaH_FILE   "DATA/Low-Res/opacCaH.dat"
#define Fe_FILE   "DATA/Low-Res/opacFe.dat"
#define K_FILE   "DATA/Low-Res/opacK.dat"
#define O3_FILE   "DATA/Low-Res/opacO3.dat"
#define Ca_FILE   "DATA/Low-Res/opacCa.dat"
#define Na_FILE   "DATA/Low-Res/opacNa.dat"
#define AlO_FILE   "DATA/Low-Res/opacAlO.dat"
#define C2H2_FILE   "DATA/Low-Res/opacC2H2.dat"
#endif /* !__INPUT_H__ */

/* ------- end ---------------------------- input.h  ----------------- */
