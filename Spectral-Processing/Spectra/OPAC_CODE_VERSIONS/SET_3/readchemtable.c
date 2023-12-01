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

  chem.el = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.el[i] = malloc(NTEMP*sizeof(double));

  chem.map_abund = malloc(NPRESSURE*sizeof(double));
  for(i=0; i<NPRESSURE; i++)
    chem.map_abund[i] = malloc(NTEMP*sizeof(double));

  /* Read in chemistry table */

  f1 = fopen(CHEM_FILE,"r");
  if(f1 == NULL){
      printf("\nreadchemtable.c:\nError opening file: No such file or directory\n\n");
      exit(1);
    }

  for (i=0; i<CHEM_FILE_NCOLS; i++) {
    fscanf(f1,"%s", dum);
  }
  // C C1H4 C1O1 C1O2 H H2 H2O1 He N N2 H3N1 O O2 O3 e-
  for (i=NPRESSURE-1; i>=0; i--)
    {
    fscanf(f1,"%le", &chem.P[i]);
    for (j=NTEMP-1; j>=0; j--)
      {
      fscanf(f1,"%le", &chem.T[j]);
      k = 1;
      k = ReadChemLine(f1, i, j, k, chem.total);
      k = ReadChemLine(f1, i, j, k, chem.C);
      k = ReadChemLine(f1, i, j, k, chem.CH4);
      k = ReadChemLine(f1, i, j, k, chem.CO);
      k = ReadChemLine(f1, i, j, k, chem.CO2);
      k = ReadChemLine(f1, i, j, k, chem.H);
      k = ReadChemLine(f1, i, j, k, chem.H2);
      k = ReadChemLine(f1, i, j, k, chem.H2O);
      k = ReadChemLine(f1, i, j, k, chem.He);
      k = ReadChemLine(f1, i, j, k, chem.N);
      k = ReadChemLine(f1, i, j, k, chem.N2);
      k = ReadChemLine(f1, i, j, k, chem.NH3);
      k = ReadChemLine(f1, i, j, k, chem.O);
      k = ReadChemLine(f1, i, j, k, chem.O2);
      k = ReadChemLine(f1, i, j, k, chem.O3);
      k = ReadChemLine(f1, i, j, k, chem.el);
      }
  }

  printf("Read in chemtable\n");
  fclose(f1);
  printf("Chemistry: \nP_0\t%e \nT_0\t%e \ntotal00 \t%e \nC\t%e \nCH4\t%e \nCO\t%e \n",
	 chem.P[0], chem.T[0], chem.total[0][0], chem.C[0][0], chem.CH4[0][0], chem.CO[0][0]);
  return;

}
/* ------- end -------------- ReadChemTable.c -------------------- */

/* ------- start ------------ FreeChemTable.c -------------------- */

/* Added: free Fe, Fe+ */
void FreeChemTable()
{
  free(chem.P);
  free(chem.T);
  free(chem.total);
  free(chem.H2);
  free(chem.H);
  free(chem.He);
  free(chem.H2O);
  free(chem.CH4);
  free(chem.CO);
  free(chem.CO2);
  free(chem.O);
  free(chem.C);
  free(chem.N);
  free(chem.NH3);
  free(chem.N2);
  free(chem.O2);
  free(chem.O3);
  free(chem.el);
  free(chem.map_abund);
}

/* ------- end -------------- FreeChemTable.c -------------------- */
