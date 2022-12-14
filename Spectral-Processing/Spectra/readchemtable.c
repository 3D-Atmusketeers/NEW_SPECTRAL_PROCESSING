
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

void ReadChemTable() {

  int i, j;
  char dum[8];

  FILE *f1;

  /* Allocate memory for Chem structure */

  chem.T = malloc(NTEMP*sizeof(double));
  chem.P = malloc(NPRESSURE*sizeof(double));

  chem.total = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.total[i] = malloc(NTEMP*sizeof(double));

  chem.H2 = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.H2[i] = malloc(NTEMP*sizeof(double));

  chem.H = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.H[i] = malloc(NTEMP*sizeof(double));

  chem.He = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.He[i] = malloc(NTEMP*sizeof(double));

  chem.H2O = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.H2O[i] = malloc(NTEMP*sizeof(double));

  chem.CH4 = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.CH4[i] = malloc(NTEMP*sizeof(double));

  chem.CO = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.CO[i] = malloc(NTEMP*sizeof(double));

  chem.CO2 = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.CO2[i] = malloc(NTEMP*sizeof(double));

  chem.O = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.O[i] = malloc(NTEMP*sizeof(double));

  chem.C = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.C[i] = malloc(NTEMP*sizeof(double));

  chem.N = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.N[i] = malloc(NTEMP*sizeof(double));

  chem.NH3 = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.NH3[i] = malloc(NTEMP*sizeof(double));

  chem.N2 = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.N2[i] = malloc(NTEMP*sizeof(double));

  chem.O2 = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.O2[i] = malloc(NTEMP*sizeof(double));

  chem.O3 = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.O3[i] = malloc(NTEMP*sizeof(double));

  /* Read in chemistry table */

  f1 = fopen(CHEM_FILE,"r");
  if(f1 == NULL){
      printf("\nreadchemtable.c:\nError opening file: No such file or directory\n\n");
      exit(1);
    }

  for (i=0; i<16; i++)
    fscanf(f1,"%s", dum);
  

  for (i=NPRESSURE-1; i>=0; i--) {
    fscanf(f1,"%le", &chem.P[i]);
    for (j=NTEMP-1; j>=0; j--) {
      fscanf(f1,"%le", &chem.T[j]);
      fscanf(f1,"%le", &chem.total[i][j]);
      fscanf(f1,"%le", &chem.C[i][j]);
      fscanf(f1,"%le", &chem.CH4[i][j]);
      fscanf(f1,"%le", &chem.CO[i][j]);
      fscanf(f1,"%le", &chem.CO2[i][j]);
      fscanf(f1,"%le", &chem.H[i][j]);
      fscanf(f1,"%le", &chem.H2[i][j]);
      fscanf(f1,"%le", &chem.H2O[i][j]);
      fscanf(f1,"%le", &chem.He[i][j]);
      fscanf(f1,"%le", &chem.N[i][j]);
      fscanf(f1,"%le", &chem.N2[i][j]);
      fscanf(f1,"%le", &chem.NH3[i][j]);
      fscanf(f1,"%le", &chem.O[i][j]);
      fscanf(f1,"%le", &chem.O2[i][j]);
      fscanf(f1,"%le", &chem.O3[i][j]);
    }
  }
  
  fclose(f1);
  printf("Chemistry: \nP_0\t%e \nT_0\t%e \ntotal\t%e \nH2\t%e \nH\t%e \nHe\t%e \nH2O\t%e \nCH4\t%e \nCO\t%e \nCO2\t%e \nO\t%e \nC\t%e \nN\t%e \nNH3\t%e \nN2\t%e \nO2\t%e \nO3\t%e \n", 
	 chem.P[0], chem.T[0], chem.total[0][0], chem.H2[0][0], chem.H[0][0], 
	 chem.He[0][0], chem.H2O[0][0], chem.CH4[0][0], chem.CO[0][0], 
	 chem.CO2[0][0], chem.O[0][0], chem.C[0][0], chem.N[0][0], 
	 chem.NH3[0][0], chem.N2[0][0], chem.O2[0][0], chem.O3[0][0]);


  return;

}

/* ------- end -------------- ReadChemTable.c -------------------- */

/* ------- start ------------ FreeChemTable.c -------------------- */

void FreeChemTable()
{
 
  free(chem.P);
  free(chem.T);
  free(chem.total);
  free(chem.H2);
  free(chem.H2O);
  free(chem.CH4);
  free(chem.CO);
  free(chem.CO2);
  free(chem.H);
  free(chem.O);
  free(chem.C);
  free(chem.N);
  free(chem.NH3);
  free(chem.N2);
  free(chem.He);
 
}

/* ------- end -------------- FreeChemTable.c -------------------- */
