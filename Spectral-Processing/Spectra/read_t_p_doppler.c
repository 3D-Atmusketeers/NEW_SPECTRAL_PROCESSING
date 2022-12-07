/*------------ file ------- readTP.c -----------------------------

Author: Eliza Miller-Ricci (emillerricci@cfa.harvard.edu)
Last modified: October 20, 2009

------------------------------------------------------------------ */

/* Reads in the temperature - pressure profile from the file 
   T_P_FILE defined in input.h
------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>

#include "input.h"
#include "atmos.h"
#include "nrutil.h"

/* --- Global variables ------------------------------------------ */

extern struct Atmos atmos;

/* --- Function prototypes --------------------------------------- */

/* ------- begin --------------------- ReadTP.c ------------------ */

void ReadTP_3D()
{
  double dum;
  double num;
  int i, j, k;
  FILE *file;

  /* Allocate Memory */

  atmos.lat = dvector(0, NLAT-1);
  atmos.lon = dvector(0, NLON-1);
  atmos.alt = dvector(0, NTAU-1);

  atmos.P_3d = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
    atmos.P_3d[i] = malloc(NLON*sizeof(double));
    for(j=0; j<NLON; j++){
      atmos.P_3d[i][j] = malloc(NTAU*sizeof(double));
    }
  }
  atmos.T_3d = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
    atmos.T_3d[i] = malloc(NLON*sizeof(double));
    for(j=0; j<NLON; j++){
      atmos.T_3d[i][j] = malloc(NTAU*sizeof(double));
    }
  }
  atmos.vel_ew = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
    atmos.vel_ew[i] = malloc(NLON*sizeof(double));
    for(j=0; j<NLON; j++){
      atmos.vel_ew[i][j] = malloc(NTAU*sizeof(double));
    }
  }
  atmos.vel_ns = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
    atmos.vel_ns[i] = malloc(NLON*sizeof(double));
    for(j=0; j<NLON; j++){
      atmos.vel_ns[i][j] = malloc(NTAU*sizeof(double));
    }
  }
  atmos.vel_ve = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
    atmos.vel_ve[i] = malloc(NLON*sizeof(double));
    for(j=0; j<NLON; j++){
      atmos.vel_ve[i][j] = malloc(NTAU*sizeof(double));
    }
  }
  /* allocate memory for aero properties */

  /* MnS */
  atmos.aero_tau_pre_qext_1 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.aero_tau_pre_qext_1[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.aero_tau_pre_qext_1[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.sw_asym_1 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_asym_1[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_asym_1[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.sw_pi0_1 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_pi0_1[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_pi0_1[i][j] = malloc(NTAU*sizeof(double));
      }
  }

  /* Al2O3 */
  atmos.aero_tau_pre_qext_2 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.aero_tau_pre_qext_2[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.aero_tau_pre_qext_2[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.sw_asym_2 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_asym_2[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_asym_2[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.sw_pi0_2 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_pi0_2[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_pi0_2[i][j] = malloc(NTAU*sizeof(double));
      }
  }

  /* Fe */
  atmos.aero_tau_pre_qext_3 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.aero_tau_pre_qext_3[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.aero_tau_pre_qext_3[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.sw_asym_3 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_asym_3[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_asym_3[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.sw_pi0_3 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_pi0_3[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_pi0_3[i][j] = malloc(NTAU*sizeof(double));
      }
  }

  /* MgSiO3 */
  atmos.aero_tau_pre_qext_4 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.aero_tau_pre_qext_4[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.aero_tau_pre_qext_4[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.sw_asym_4 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_asym_4[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_asym_4[i][j] = malloc(NTAU*sizeof(double));
      }
  }

  atmos.sw_pi0_4 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_pi0_4[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_pi0_4[i][j] = malloc(NTAU*sizeof(double));
      }
  }




  atmos.aero_tau_pre_qext_5 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.aero_tau_pre_qext_5[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.aero_tau_pre_qext_5[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.sw_asym_5 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_asym_5[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_asym_5[i][j] = malloc(NTAU*sizeof(double));
      }
  }

  atmos.sw_pi0_5 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_pi0_5[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_pi0_5[i][j] = malloc(NTAU*sizeof(double));
      }
  }

  atmos.aero_tau_pre_qext_6 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.aero_tau_pre_qext_6[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.aero_tau_pre_qext_6[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.sw_asym_6 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_asym_6[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_asym_6[i][j] = malloc(NTAU*sizeof(double));
      }
  }

  atmos.sw_pi0_6 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_pi0_6[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_pi0_6[i][j] = malloc(NTAU*sizeof(double));
      }
  }

  atmos.aero_tau_pre_qext_7 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.aero_tau_pre_qext_7[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.aero_tau_pre_qext_7[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.sw_asym_7 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_asym_7[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_asym_7[i][j] = malloc(NTAU*sizeof(double));
      }
  }

  atmos.sw_pi0_7 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_pi0_7[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_pi0_7[i][j] = malloc(NTAU*sizeof(double));
      }
  }


  atmos.aero_tau_pre_qext_8 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.aero_tau_pre_qext_8[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.aero_tau_pre_qext_8[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.sw_asym_8 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_asym_8[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_asym_8[i][j] = malloc(NTAU*sizeof(double));
      }
  }

  atmos.sw_pi0_8 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_pi0_8[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_pi0_8[i][j] = malloc(NTAU*sizeof(double));
      }
  }



  atmos.aero_tau_pre_qext_9 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.aero_tau_pre_qext_9[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.aero_tau_pre_qext_9[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.sw_asym_9 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_asym_9[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_asym_9[i][j] = malloc(NTAU*sizeof(double));
      }
  }

  atmos.sw_pi0_9 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_pi0_9[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_pi0_9[i][j] = malloc(NTAU*sizeof(double));
      }
  }



  atmos.aero_tau_pre_qext_10 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.aero_tau_pre_qext_10[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.aero_tau_pre_qext_10[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.sw_asym_10 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_asym_10[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_asym_10[i][j] = malloc(NTAU*sizeof(double));
      }
  }

  atmos.sw_pi0_10 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_pi0_10[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_pi0_10[i][j] = malloc(NTAU*sizeof(double));
      }
  }

  atmos.aero_tau_pre_qext_11 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.aero_tau_pre_qext_11[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.aero_tau_pre_qext_11[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.sw_asym_11 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_asym_11[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_asym_11[i][j] = malloc(NTAU*sizeof(double));
      }
  }

  atmos.sw_pi0_11 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_pi0_11[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_pi0_11[i][j] = malloc(NTAU*sizeof(double));
      }
  }

  atmos.aero_tau_pre_qext_12 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.aero_tau_pre_qext_12[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.aero_tau_pre_qext_12[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.sw_asym_12 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_asym_12[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_asym_12[i][j] = malloc(NTAU*sizeof(double));
      }
  }

  atmos.sw_pi0_12 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_pi0_12[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_pi0_12[i][j] = malloc(NTAU*sizeof(double));
      }
  }

  atmos.aero_tau_pre_qext_13 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.aero_tau_pre_qext_13[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.aero_tau_pre_qext_13[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.sw_asym_13 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_asym_13[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_asym_13[i][j] = malloc(NTAU*sizeof(double));
      }
  }

  atmos.sw_pi0_13 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.sw_pi0_13[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.sw_pi0_13[i][j] = malloc(NTAU*sizeof(double));
      }
  }


  atmos.aero_tau_haze = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.aero_tau_haze[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.aero_tau_haze[i][j] = malloc(NTAU*sizeof(double));
      }
  }
  atmos.tau_asym = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.tau_asym[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.tau_asym[i][j] = malloc(NTAU*sizeof(double));
      }
  }

  atmos.tau_pi0 = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
      atmos.tau_pi0[i] = malloc(NLON*sizeof(double));
      for(j=0; j<NLON; j++){
          atmos.tau_pi0[i][j] = malloc(NTAU*sizeof(double));
      }
  }


  atmos.incident_frac = malloc(NLAT*sizeof(double));
  for(i=0; i<NLAT; i++){
    atmos.incident_frac[i] = malloc(NLON*sizeof(double));
    for(j=0; j<NLON; j++){
      atmos.incident_frac[i][j] = malloc(NTAU*sizeof(double));
    }
  }






    file = fopen(T_P_3D_FILE, "r");
    if(file == NULL){
        printf("\nreadt_p.c:\nError opening file: No such file or directory\n\n");
        exit(1);
    }

    for (i=0; i<NLAT; i++){
        for(j=0; j<NLON; j++){
            for(k=0; k<NTAU; k++)
            {
                if (fscanf(file, "%le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le",
                       &atmos.lat[i], &atmos.lon[j], &num,
                       &atmos.alt[k], &atmos.P_3d[i][j][k], &atmos.T_3d[i][j][k],
                       &atmos.vel_ew[i][j][k], &atmos.vel_ns[i][j][k], &atmos.vel_ve[i][j][k],
                       &atmos.aero_tau_pre_qext_1[i][j][k], &atmos.sw_asym_1[i][j][k], &atmos.sw_pi0_1[i][j][k],
                       &atmos.aero_tau_pre_qext_2[i][j][k], &atmos.sw_asym_2[i][j][k], &atmos.sw_pi0_2[i][j][k],
                       &atmos.aero_tau_pre_qext_3[i][j][k], &atmos.sw_asym_3[i][j][k], &atmos.sw_pi0_3[i][j][k],
                       &atmos.aero_tau_pre_qext_4[i][j][k], &atmos.sw_asym_4[i][j][k], &atmos.sw_pi0_4[i][j][k],
                       &atmos.aero_tau_pre_qext_5[i][j][k], &atmos.sw_asym_5[i][j][k], &atmos.sw_pi0_5[i][j][k],
                       &atmos.aero_tau_pre_qext_6[i][j][k], &atmos.sw_asym_6[i][j][k], &atmos.sw_pi0_6[i][j][k],
                       &atmos.aero_tau_pre_qext_7[i][j][k], &atmos.sw_asym_7[i][j][k], &atmos.sw_pi0_7[i][j][k],
                       &atmos.aero_tau_pre_qext_8[i][j][k], &atmos.sw_asym_8[i][j][k], &atmos.sw_pi0_8[i][j][k],
                       &atmos.aero_tau_pre_qext_9[i][j][k], &atmos.sw_asym_9[i][j][k], &atmos.sw_pi0_9[i][j][k],
                       &atmos.aero_tau_pre_qext_10[i][j][k], &atmos.sw_asym_10[i][j][k], &atmos.sw_pi0_10[i][j][k],
                       &atmos.aero_tau_pre_qext_11[i][j][k], &atmos.sw_asym_11[i][j][k], &atmos.sw_pi0_11[i][j][k],
                       &atmos.aero_tau_pre_qext_12[i][j][k], &atmos.sw_asym_12[i][j][k], &atmos.sw_pi0_12[i][j][k],
                       &atmos.aero_tau_pre_qext_13[i][j][k], &atmos.sw_asym_13[i][j][k], &atmos.sw_pi0_13[i][j][k],
                       &atmos.aero_tau_haze[i][j][k], &atmos.tau_asym[i][j][k], &atmos.tau_pi0[i][j][k],
                       &atmos.incident_frac[i][j][k]));
            }



            for(k=NTAU; k<NTAU; k++)
            {
                if (fscanf(file, "%le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le %le %le \
                                  %le",
                                  &dum, &dum, &num,
                                  &dum, &dum, &dum,
                                  &dum, &dum, &dum,
                                  &dum, &dum, &dum,
                                  &dum, &dum, &dum,
                                  &dum, &dum, &dum,
                                  &dum, &dum, &dum,
                                  &dum, &dum, &dum,
                                  &dum, &dum, &dum,
                                  &dum, &dum, &dum,
                                  &dum, &dum, &dum,
                                  &dum, &dum, &dum,
                                  &dum, &dum, &dum,
                                  &dum, &dum, &dum,
                                  &dum, &dum, &dum,
                                  &dum, &dum, &dum,
                                  &dum, &dum, &dum,
                                  &dum));
            }

        }
    }

    printf("The First Profile, 5 layers\n");
    for(k=0; k<5; k++)
    {
        printf("%.2e %.2e %d %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e\n",
                atmos.lat[0], atmos.lon[0], k + 1, \
                atmos.alt[k], atmos.P_3d[0][0][k], atmos.T_3d[0][0][k], \
                atmos.vel_ew[0][0][k], atmos.vel_ns[0][0][k], atmos.vel_ve[0][0][k], \
                atmos.aero_tau_pre_qext_4[0][0][k], atmos.sw_asym_4[0][0][k], atmos.sw_pi0_4[0][0][k], \
                atmos.aero_tau_pre_qext_5[0][0][k], atmos.sw_asym_5[0][0][k], atmos.sw_pi0_5[0][0][k], \
                atmos.aero_tau_haze[0][0][k], atmos.tau_asym[0][0][k], atmos.tau_pi0[0][0][k], \
                atmos.incident_frac[0][0][k]);
    }


            //        %.1e %.1e %.1e \
        //        %.1e %.1e %.1e \
        //        %.1e %.1e %.1e \
        //        %.1e %.1e %.1e \
        //        %.1e %.1e %.1e \
        //        %.1e %.1e %.1e \
        //        %.1e %.1e %.1e \
        //        %.1e %.1e %.1e \
        //        %.1e %.1e %.1e \
        //        %.1e %.1e %.1e \
        //        %.1e %.1e %.1e \
        //        %.1e %.1e %.1e \

            //            atmos.aero_tau_pre_qext_3[0][0][k], atmos.sw_asym_3[0][0][k], atmos.sw_pi0_3[0][0][k],\
    //            atmos.aero_tau_pre_qext_4[0][0][k], atmos.sw_asym_4[0][0][k], atmos.sw_pi0_4[0][0][k],\
    //            atmos.aero_tau_pre_qext_5[0][0][k], atmos.sw_asym_5[0][0][k], atmos.sw_pi0_5[0][0][k],\
    //            atmos.aero_tau_pre_qext_6[0][0][k], atmos.sw_asym_6[0][0][k], atmos.sw_pi0_6[0][0][k],\
    //            atmos.aero_tau_pre_qext_7[0][0][k], atmos.sw_asym_7[0][0][k], atmos.sw_pi0_7[0][0][k],\
    //            atmos.aero_tau_pre_qext_8[0][0][k], atmos.sw_asym_8[0][0][k], atmos.sw_pi0_8[0][0][k],\
    //            atmos.aero_tau_pre_qext_9[0][0][k], atmos.sw_asym_9[0][0][k], atmos.sw_pi0_9[0][0][k],\
    //            atmos.aero_tau_pre_qext_10[0][0][k], atmos.sw_asym_10[0][0][k], atmos.sw_pi0_10[0][0][k],\
    //            atmos.aero_tau_pre_qext_11[0][0][k], atmos.sw_asym_11[0][0][k], atmos.sw_pi0_11[0][0][k],\
    //            atmos.aero_tau_pre_qext_12[0][0][k], atmos.sw_asym_12[0][0][k], atmos.sw_pi0_12[0][0][k],\
    //           atmos.aero_tau_pre_qext_13[0][0][k], atmos.sw_asym_13[0][0][k], atmos.sw_pi0_13[0][0][k],\

    fclose(file);

}

/* ------- end ----------------------- ReadTP.c ------------------ */
