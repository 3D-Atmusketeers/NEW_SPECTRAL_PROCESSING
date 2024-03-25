/*----------------------- readchemtable.c ------------------------
Author: Eliza Miller-Ricci (emillerricci@cfa.harvard.edu)
Last modified: October 20, 2009
------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>

#include "input.h"
#include "opac.h"

/* --- Global variables ------------------------------------------ */

extern struct Chem chem;

/* --- Function prototypes --------------------------------------- */

/* ---------------------------------------------------------------
 * Read in chemistry files: abundance(pressure, temperature)
 * --------------------------------------------------------------- */

/* ------- begin ------------ ReadChemTable.c -------------------- */


int ReadChemLine(FILE *f1, int i, int j, int k, double **chem_array) {
    /* Reads a chemistry line. keeps a running tally k of how many
    columns have been read in already. Also takes care of reading into
    the map_abund array. */
    fscanf(f1,"%le", &chem_array[i][j]);
    if (ABUND_TRACK_IND <= 1) {
        printf("\nMapping temperature/pressure doesn't make sense here!\n");
        exit(1);
    }
    if (k == ABUND_TRACK_IND) {
        chem.map_abund[i][j] = chem_array[i][j];
    }

    k = k + 1;
    return k;
}

void ReadChemTable()
  {
  int i, j, k;
  char dum[1000];

  FILE *f1;
  /* Allocate memory for Chem structure */

  chem.T = malloc(NTEMP*sizeof(double));
  chem.P = malloc(NPRESSURE*sizeof(double));





chem.total = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.total[i] = malloc(NTEMP*sizeof(double));


  chem.CH4 = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.CH4[i] = malloc(NTEMP*sizeof(double));


  chem.PH3 = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.PH3[i] = malloc(NTEMP*sizeof(double));


  chem.K = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.K[i] = malloc(NTEMP*sizeof(double));


  chem.CrH = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.CrH[i] = malloc(NTEMP*sizeof(double));


  chem.He = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.He[i] = malloc(NTEMP*sizeof(double));


  chem.H2S = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.H2S[i] = malloc(NTEMP*sizeof(double));


  chem.Ca_plus = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.Ca_plus[i] = malloc(NTEMP*sizeof(double));


  chem.VO = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.VO[i] = malloc(NTEMP*sizeof(double));


  chem.SiH4 = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.SiH4[i] = malloc(NTEMP*sizeof(double));


  chem.N = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.N[i] = malloc(NTEMP*sizeof(double));


  chem.Fe_plus = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.Fe_plus[i] = malloc(NTEMP*sizeof(double));


  chem.H2 = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.H2[i] = malloc(NTEMP*sizeof(double));


  chem.TiO = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.TiO[i] = malloc(NTEMP*sizeof(double));


  chem.Fe = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.Fe[i] = malloc(NTEMP*sizeof(double));


  chem.el = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.el[i] = malloc(NTEMP*sizeof(double));


  chem.CO2 = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.CO2[i] = malloc(NTEMP*sizeof(double));


  chem.CaH = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.CaH[i] = malloc(NTEMP*sizeof(double));


  chem.SiH = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.SiH[i] = malloc(NTEMP*sizeof(double));


  chem.SiO2 = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.SiO2[i] = malloc(NTEMP*sizeof(double));


  chem.C = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.C[i] = malloc(NTEMP*sizeof(double));


  chem.H2O = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.H2O[i] = malloc(NTEMP*sizeof(double));


  chem.SH = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.SH[i] = malloc(NTEMP*sizeof(double));


  chem.SiO = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.SiO[i] = malloc(NTEMP*sizeof(double));


  chem.CaO = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.CaO[i] = malloc(NTEMP*sizeof(double));


  chem.AlH = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.AlH[i] = malloc(NTEMP*sizeof(double));


  chem.Na = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.Na[i] = malloc(NTEMP*sizeof(double));


  chem.O3 = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.O3[i] = malloc(NTEMP*sizeof(double));


  chem.AlO = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.AlO[i] = malloc(NTEMP*sizeof(double));


  chem.CO = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.CO[i] = malloc(NTEMP*sizeof(double));


  chem.C2H2 = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.C2H2[i] = malloc(NTEMP*sizeof(double));


  chem.N2 = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.N2[i] = malloc(NTEMP*sizeof(double));


  chem.O2 = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.O2[i] = malloc(NTEMP*sizeof(double));


  chem.NaH = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.NaH[i] = malloc(NTEMP*sizeof(double));


  chem.Ca = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.Ca[i] = malloc(NTEMP*sizeof(double));


  chem.NH3 = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.NH3[i] = malloc(NTEMP*sizeof(double));


  chem.MgH = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.MgH[i] = malloc(NTEMP*sizeof(double));


  chem.FeH = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.FeH[i] = malloc(NTEMP*sizeof(double));


  chem.NO = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.NO[i] = malloc(NTEMP*sizeof(double));


  chem.SO2 = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.SO2[i] = malloc(NTEMP*sizeof(double));


  chem.HF = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.HF[i] = malloc(NTEMP*sizeof(double));


  chem.O = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.O[i] = malloc(NTEMP*sizeof(double));


  chem.OH = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.OH[i] = malloc(NTEMP*sizeof(double));


  chem.Mg = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.Mg[i] = malloc(NTEMP*sizeof(double));


  chem.H = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.H[i] = malloc(NTEMP*sizeof(double));


  chem.HCN = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.HCN[i] = malloc(NTEMP*sizeof(double));


  chem.map_abund = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.map_abund[i] = malloc(NTEMP*sizeof(double));


  chem.rho_k_T = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.rho_k_T[i] = malloc(NTEMP*sizeof(double));


  chem.mu = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.mu[i] = malloc(NTEMP*sizeof(double));




  /* Read in chemistry table */

  f1 = fopen(CHEM_FILE,"r");
  if(f1 == NULL){
      printf("\nreadchemtable.c:\nError opening file: No such file or directory\n\n");
      exit(1);
    }

  for (i=0; i<CHEM_FILE_NCOLS; i++) {
    fscanf(f1,"%s", dum);
  }
  for (i=NPRESSURE-1; i>=0; i--)
    {
    fscanf(f1,"%le", &chem.P[i]);
    for (j=NTEMP-1; j>=0; j--)
      {
      fscanf(f1,"%le", &chem.T[j]);
      k = 1;
      k = ReadChemLine(f1, i, j, k, chem.total);

      k = ReadChemLine(f1, i, j, k, chem.CH4);
      k = ReadChemLine(f1, i, j, k, chem.PH3);
      k = ReadChemLine(f1, i, j, k, chem.K);
      k = ReadChemLine(f1, i, j, k, chem.CrH);
      k = ReadChemLine(f1, i, j, k, chem.He);
      k = ReadChemLine(f1, i, j, k, chem.H2S);
      k = ReadChemLine(f1, i, j, k, chem.Ca_plus);
      k = ReadChemLine(f1, i, j, k, chem.VO);
      k = ReadChemLine(f1, i, j, k, chem.SiH4);
      k = ReadChemLine(f1, i, j, k, chem.N);
      k = ReadChemLine(f1, i, j, k, chem.Fe_plus);
      k = ReadChemLine(f1, i, j, k, chem.H2);
      k = ReadChemLine(f1, i, j, k, chem.TiO);
      k = ReadChemLine(f1, i, j, k, chem.Fe);
      k = ReadChemLine(f1, i, j, k, chem.el);
      k = ReadChemLine(f1, i, j, k, chem.CO2);
      k = ReadChemLine(f1, i, j, k, chem.CaH);
      k = ReadChemLine(f1, i, j, k, chem.SiH);
      k = ReadChemLine(f1, i, j, k, chem.SiO2);
      k = ReadChemLine(f1, i, j, k, chem.C);
      k = ReadChemLine(f1, i, j, k, chem.H2O);
      k = ReadChemLine(f1, i, j, k, chem.SH);
      k = ReadChemLine(f1, i, j, k, chem.SiO);
      k = ReadChemLine(f1, i, j, k, chem.CaO);
      k = ReadChemLine(f1, i, j, k, chem.AlH);
      k = ReadChemLine(f1, i, j, k, chem.Na);
      k = ReadChemLine(f1, i, j, k, chem.O3);
      k = ReadChemLine(f1, i, j, k, chem.AlO);
      k = ReadChemLine(f1, i, j, k, chem.CO);
      k = ReadChemLine(f1, i, j, k, chem.C2H2);
      k = ReadChemLine(f1, i, j, k, chem.N2);
      k = ReadChemLine(f1, i, j, k, chem.O2);
      k = ReadChemLine(f1, i, j, k, chem.NaH);
      k = ReadChemLine(f1, i, j, k, chem.Ca);
      k = ReadChemLine(f1, i, j, k, chem.NH3);
      k = ReadChemLine(f1, i, j, k, chem.MgH);
      k = ReadChemLine(f1, i, j, k, chem.FeH);
      k = ReadChemLine(f1, i, j, k, chem.NO);
      k = ReadChemLine(f1, i, j, k, chem.SO2);
      k = ReadChemLine(f1, i, j, k, chem.HF);
      k = ReadChemLine(f1, i, j, k, chem.O);
      k = ReadChemLine(f1, i, j, k, chem.OH);
      k = ReadChemLine(f1, i, j, k, chem.Mg);
      k = ReadChemLine(f1, i, j, k, chem.H);
      k = ReadChemLine(f1, i, j, k, chem.HCN);

      //k = ReadChemLine(f1, i, j, k, chem.el);
      //k = ReadChemLine(f1, i, j, k, chem.H);
      //k = ReadChemLine(f1, i, j, k, chem.H2);
      //k = ReadChemLine(f1, i, j, k, chem.He);
      //k = ReadChemLine(f1, i, j, k, chem.C2H2);
      //k = ReadChemLine(f1, i, j, k, chem.CH4);
      //k = ReadChemLine(f1, i, j, k, chem.CO);
      //k = ReadChemLine(f1, i, j, k, chem.CO2);
      //k = ReadChemLine(f1, i, j, k, chem.H2O);
      //k = ReadChemLine(f1, i, j, k, chem.H2S);
      //k = ReadChemLine(f1, i, j, k, chem.HCN);
      //k = ReadChemLine(f1, i, j, k, chem.K);
      //k = ReadChemLine(f1, i, j, k, chem.Na);
      //k = ReadChemLine(f1, i, j, k, chem.NH3);
      //k = ReadChemLine(f1, i, j, k, chem.PH3);
      //k = ReadChemLine(f1, i, j, k, chem.TiO);
      //k = ReadChemLine(f1, i, j, k, chem.VO);
      //k = ReadChemLine(f1, i, j, k, chem.FeH);
      //k = ReadChemLine(f1, i, j, k, chem.Ar);
      //k = ReadChemLine(f1, i, j, k, chem.N2);
      //k = ReadChemLine(f1, i, j, k, chem.O2);
      }
  }
  fclose(f1);

  printf("======================================================");
  printf("\n==== Debug Chem Information ====\n");
  printf("Read in chemtable.\n");
  printf("Printing out an example:\n");
  printf("\tTemperature Range: %.2f to %.2f\n", chem.T[0], chem.T[NTEMP-1]);
  printf("\tPressure Range: %.2f to %.2f\n", chem.P[0], chem.P[NPRESSURE-1]);
  printf("\tCH4 abundance at a pressure of %.2f and a temperature of %.2f : %.4e\n", chem.P[0], chem.T[0], chem.CH4[0][0]);

  printf("======================================================\n\n");

  return;

}
/* ------- end -------------- ReadChemTable.c -------------------- */

/* ------- start ------------ FreeChemTable.c -------------------- */

void FreeChemTable()
{
  free(chem.P);
  free(chem.T);
  free(chem.total);
  free(chem.CH4);
  free(chem.PH3);
  free(chem.K);
  free(chem.CrH);
  free(chem.He);
  free(chem.H2S);
  free(chem.Ca_plus);
  free(chem.VO);
  free(chem.SiH4);
  free(chem.N);
  free(chem.Fe_plus);
  free(chem.H2);
  free(chem.TiO);
  free(chem.Fe);
  free(chem.el);
  free(chem.CO2);
  free(chem.CaH);
  free(chem.SiH);
  free(chem.SiO2);
  free(chem.C);
  free(chem.H2O);
  free(chem.SH);
  free(chem.SiO);
  free(chem.CaO);
  free(chem.AlH);
  free(chem.Na);
  free(chem.O3);
  free(chem.AlO);
  free(chem.CO);
  free(chem.C2H2);
  free(chem.N2);
  free(chem.O2);
  free(chem.NaH);
  free(chem.Ca);
  free(chem.NH3);
  free(chem.MgH);
  free(chem.FeH);
  free(chem.NO);
  free(chem.SO2);
  free(chem.HF);
  free(chem.O);
  free(chem.OH);
  free(chem.Mg);
  free(chem.H);
  free(chem.HCN);
  free(chem.map_abund);
  free(chem.rho_k_T);
  free(chem.mu);
}
/* ------- end -------------- FreeChemTable.c -------------------- */
