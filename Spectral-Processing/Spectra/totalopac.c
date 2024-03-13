/*----------------------- totalopac.c ----------------------------

Author: Eliza Miller-Ricci (emillerricci@cfa.harvard.edu)
Last modified: June 13, 2007

------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "input.h"
#include "include.h"
#include "constant.h"
#include "atmos.h"
#include "opac.h"
#include "nrutil.h"

/* --- Global variables ------------------------------------------ */

extern struct Opac opac;
extern struct Atmos atmos;

struct Chem chem;

// Add in the new lines for each opacity species
struct Opac opacH2O;

struct Opac opacscat;
struct Opac opacCIA;
 
/* --- Function prototypes --------------------------------------- */

void Locate(int n, double *array, double value, int *ilow);
void ReadOpacTable(struct Opac opac, char *filename);
double lint2D(double x1, double x2, double y1, double y2, double z1,
              double z2, double z3, double z4, double x, double y);
void FreeOpacTable(struct Opac opac);
void ReadChemTable();
void FreeChemTable();

/* ---------------------------------------------------------------
 * Computes the total opacity due to all of the atmospheric 
 * constituents.
 * --------------------------------------------------------------- */

/* ------- begin ------------ TotalOpac.c ------------------------ */

void TotalOpac() {
  double **opac_CIA_H2H2, **opac_CIA_H2He, **opac_CIA_H2H, **opac_CIA_H2CH4, **opac_CIA_CH4Ar,
         **opac_CIA_CH4CH4, **opac_CIA_CO2CO2, **opac_CIA_HeH, **opac_CIA_N2CH4, **opac_CIA_N2H2,
         **opac_CIA_N2N2, **opac_CIA_O2CO2, **opac_CIA_O2N2, **opac_CIA_O2O2, **opac_CIA_Hel;
  double *t_CIA, *lambda_CIA, **opac_CIA;
  int i, j, k, a, b, dum;
  char filename[65];

  FILE *f1;

  /* Allocate Memory */
  opac_CIA_H2H2 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_H2He = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_H2H = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_H2CH4 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_CH4Ar = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_CH4CH4 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_CO2CO2 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_HeH = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_N2CH4 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_N2H2 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_N2N2 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_O2CO2 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_O2N2 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_O2O2 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_Hel = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);


  t_CIA = dvector(0, 18);
  lambda_CIA = dvector(0, 999);
  opac_CIA = dmatrix(0, 19, 0, 999);

  /* Read Chemistry Table */
  ReadChemTable();
  printf("ReadChemTable done\n");

  /* Fill in mean molecular weight (mu) values */
  opac.mu = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for (j=0; j<NPRESSURE; j++) {
    for (k=0; k<NTEMP; k++)
    {
      opac.mu[j][k] = MU;
    }
  }

  /* Fill in scattering coefficients */
  opacscat.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacscat.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacscat.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }



  // Fill in the opacities for each species
  /* Fill in H2O opacities */
  opacH2O.T = dvector(0, NTEMP-1);
  opacH2O.P = dvector(0, NPRESSURE-1);
  opacH2O.Plog10 = dvector(0, NPRESSURE-1);
  opacH2O.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacH2O.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacH2O.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacH2O.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacH2O.abundance[j][k] = chem.H2O[j][k];
    }
  }
  strcpy(filename, H2O_FILE);
  ReadOpacTable(opacH2O, filename);
  printf("Read H2O Opacity done\n");




  /* Fill in collision-induced opacities */
  //f1 = fopen(CIA_FILE, "r");
  f1 = fopen("DATA/Low-Res/opacCIA.dat", "r");

  if(f1 == NULL){
    printf("\n totalopac.c:\nError opening CIA file: -- No such file or directory\n\n");
  }
  printf("CIA file found, good job!\n");
  opacCIA.T = dvector(0, NTEMP-1);
  opacCIA.P = dvector(0, NPRESSURE-1);
  opacCIA.Plog10 = dvector(0, NPRESSURE-1);



  for(i=0; i<NTEMP; i++){
    fscanf(f1, "%le", &opacCIA.T[i]);
  }
  for (k=0; k<NTEMP; k++){
    fscanf(f1, "%le", &opacCIA.T[k]);
    for (i=0; i<NLAMBDA; i++){
      fscanf(f1, "%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le", &atmos.lambda[i],
                                                                                    &opac_CIA_Hel[k][i],
                                                                                    &opac_CIA_H2H2[k][i],
                                                                                    &opac_CIA_H2He[k][i],
                                                                                    &opac_CIA_H2H[k][i],
                                                                                    &opac_CIA_H2CH4[k][i],
                                                                                    &opac_CIA_CH4Ar[k][i],
                                                                                    &opac_CIA_CH4CH4[k][i],
                                                                                    &opac_CIA_CO2CO2[k][i],
                                                                                    &opac_CIA_HeH[k][i],
                                                                                    &opac_CIA_N2CH4[k][i],
                                                                                    &opac_CIA_N2H2[k][i],
                                                                                    &opac_CIA_N2N2[k][i],
                                                                                    &opac_CIA_O2CO2[k][i],
                                                                                    &opac_CIA_O2N2[k][i],
                                                                                    &opac_CIA_O2O2[k][i]);
    }
  }


  printf("\n==== Debug Information for Lambda and Hel Columns ====\n");
  printf("Lambda Range: %e to %e\n", atmos.lambda[0], atmos.lambda[NLAMBDA-1]);
  printf("First 4 Lambdas: %e, %e, %e, %e\n", atmos.lambda[0], atmos.lambda[1], atmos.lambda[2], atmos.lambda[3]);
  printf("Hel Opacity for First Temperature at First 4 Lambdas:\n");
  printf("\tLambda 1: %0.7e\n", opac_CIA_Hel[0][0]);
  printf("\tLambda 2: %0.7e\n", opac_CIA_Hel[0][1]);
  printf("\tLambda 3: %0.7e\n", opac_CIA_Hel[0][2]);
  printf("\tLambda 4: %0.7e\n", opac_CIA_Hel[0][3]);
  printf("\tLambda 5: %0.7e\n", opac_CIA_Hel[0][4]);
  printf("=====================================================\n\n");


  printf("Stopping here\n");
  exit(0);

  // Allocate memory for opacCIA.kappa
  opacCIA.kappa = (double ***)malloc(NLAMBDA * sizeof(double **));
  if (opacCIA.kappa == NULL) {
      printf("Failed to allocate memory for opacCIA.kappa (outer dimension)\n");
      // Handle the error appropriately, e.g., return or exit the function
  }
  for (i = 0; i < NLAMBDA; i++) {
      opacCIA.kappa[i] = (double **)malloc(NPRESSURE * sizeof(double *));
      if (opacCIA.kappa[i] == NULL) {
          printf("Failed to allocate memory for opacCIA.kappa[%d] (middle dimension)\n", i);
          // Handle the error appropriately, e.g., free previously allocated memory and return or exit the function
      }
      for (j = 0; j < NPRESSURE; j++) {
          opacCIA.kappa[i][j] = (double *)malloc(NTEMP * sizeof(double));
          if (opacCIA.kappa[i][j] == NULL) {
              printf("Failed to allocate memory for opacCIA.kappa[%d][%d] (inner dimension)\n", i, j);
              // Handle the error appropriately, e.g., free previously allocated memory and return or exit the function
          }
      }
  }

  for (i=0; i<NLAMBDA; i++){
    for (j=0; j<NPRESSURE; j++){
      for (k=0; k<NTEMP; k++){

      opacCIA.kappa[i][j][k] += opac_CIA_H2H2[k][i] *
        (chem.H2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k])) *
        (chem.H2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]));

      opacCIA.kappa[i][j][k] += opac_CIA_H2He[k][i] *
        chem.H2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]) *
        chem.He[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]);

      opacCIA.kappa[i][j][k] += opac_CIA_H2H[k][i] *
        chem.H2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]) *
        chem.H[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]);

      opacCIA.kappa[i][j][k] += opac_CIA_H2CH4[k][i] *
        chem.H2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]) *
        chem.CH4[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]);

      opacCIA.kappa[i][j][k] += opac_CIA_CH4CH4[k][i] *
        chem.CH4[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]) *
        chem.CH4[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]);

      opacCIA.kappa[i][j][k] += opac_CIA_CO2CO2[k][i] *
        chem.CO2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]) *
        chem.CO2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]);

      opacCIA.kappa[i][j][k] += opac_CIA_HeH[k][i] *
        chem.He[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]) *
        chem.H[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]);

      opacCIA.kappa[i][j][k] += opac_CIA_N2CH4[k][i] *
        chem.N2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]) *
        chem.CH4[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]);

      opacCIA.kappa[i][j][k] += opac_CIA_N2H2[k][i] *
        chem.N2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]) *
        chem.H2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]);

      opacCIA.kappa[i][j][k] += opac_CIA_N2N2[k][i] *
        chem.N2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]) *
        chem.N2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]);

      opacCIA.kappa[i][j][k] += opac_CIA_O2CO2[k][i] *
        chem.O2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]) *
        chem.CO2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]);

      opacCIA.kappa[i][j][k] += opac_CIA_O2N2[k][i] *
        chem.O2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]) *
        chem.N2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]);

      opacCIA.kappa[i][j][k] += opac_CIA_O2O2[k][i] *
        chem.O2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]) *
        chem.O2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]);

      opacCIA.kappa[i][j][k] += opac_CIA_Hel[k][i] *
	       chem.H[j][k] * chem.P[j] / (KBOLTZMANN * chem.T[k]) *  
	       chem.el[j][k] * chem.P[j] / (KBOLTZMANN * chem.T[k]);
      }
    }
  }
  fclose(f1);

  /* Rayleigh scattering */
  /* (Polarizabilities from the CRC Handbook) */

  for (i=0; i<NLAMBDA; i++)
    {
    for (j=0; j<NPRESSURE; j++)
      {
      for (k=0; k<NTEMP; k++)
      {
	  opacscat.kappa[i][j][k] +=
	  (8.0*PI/3.0) * SQ(0.80e-30) *
	  SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	  chem.H2[j][k]*chem.P[j] / (KBOLTZMANN * chem.T[k])
	  +
	  (8.0*PI/3.0) * SQ(0.21e-30) *
	  SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	  chem.He[j][k]*chem.P[j] / (KBOLTZMANN * chem.T[k])
	  +
	  (8.0*PI/3.0) * SQ(1.45e-30) *
	  SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	  chem.H2O[j][k]*chem.P[j] / (KBOLTZMANN * chem.T[k])
	  +
	  (8.0*PI/3.0) * SQ(1.95e-30) *
	  SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	  chem.CO[j][k]*chem.P[j] / (KBOLTZMANN * chem.T[k])
	  +
	  (8.0*PI/3.0) * SQ(2.91e-30) *
	  SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	  chem.CO2[j][k]*chem.P[j] / (KBOLTZMANN * chem.T[k])
	  +
	  (8.0*PI/3.0) * SQ(2.26e-30) *
	  SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	  chem.NH3[j][k]*chem.P[j] / (KBOLTZMANN * chem.T[k])
	  +
	  (8.0*PI/3.0) * SQ(2.59e-30) *
	  SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	  chem.CH4[j][k]*chem.P[j] / (KBOLTZMANN * chem.T[k]);
      }
    }
  }

  /* Fill in total opacities */
  opac.T = malloc(NTEMP*sizeof(double));
  opac.P = malloc(NPRESSURE*sizeof(double));
  opac.Plog10 = malloc(NPRESSURE*sizeof(double));
  opac.kappa = malloc(NLAMBDA*sizeof(double));
  opac.alpha = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opac.kappa[i] = malloc(NPRESSURE*sizeof(double));
    opac.alpha[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opac.kappa[i][j] = malloc(NTEMP*sizeof(double));
      opac.alpha[i][j] = malloc(NTEMP*sizeof(double));
    }
  }

  for (k=0; k<NTEMP; k++)
    opac.T[k] = opacH2O.T[k];

  for (j=0; j<NPRESSURE; j++) {
    opac.P[j] = opacH2O.P[j];
    opac.Plog10[j] = opacH2O.Plog10[j];
  }

  // Do a final summation over the species
  for (i=0; i<NLAMBDA; i++) {
    for (j=0; j<NPRESSURE; j++) {
      for (k=0; k<NTEMP; k++) {
          opac.kappa[i][j][k] =                              + opacH2O.kappa[i][j][k]
                              + opacscat.kappa[i][j][k]
                              + opacCIA.kappa[i][j][k];
      }
    }
  }


  
  // Free not needed opacity structures and chemistry table
  FreeOpacTable(opacH2O);


  FreeOpacTable(opacCIA);
  FreeOpacTable(opacscat);
  FreeChemTable();

  free_dmatrix(opac_CIA_H2H2, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_H2He, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_H2H, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_H2CH4, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_CH4Ar, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_CH4CH4, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_CO2CO2, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_HeH, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_N2CH4, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_N2H2, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_N2N2, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_O2CO2, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_O2N2, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_O2O2, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_Hel, 0, NTEMP-1, 0, NLAMBDA-1);

}

/* ------- end -------------- TotalOpac.c ------------------------ */
