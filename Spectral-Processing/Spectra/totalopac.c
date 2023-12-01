/*----------------------- totalopac.c ----------------------------

Author: Eliza Miller-Ricci (emillerricci@cfa.harvard.edu)

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

struct Opac opacCH4;
struct Opac opacCO2;
struct Opac opacCO;
struct Opac opacH2O;
struct Opac opacNH3;
struct Opac opacO2;
struct Opac opacO3;
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

  double **opac_CIA_H_el, **opac_CIA_He_H, **opac_CIA_CH4_CH4, **opac_CIA_H2_He, **opac_CIA_H2_CH4,
         **opac_CIA_H2_H, **opac_CIA_H2_H2, **opac_CIA_CO2_CO2;
  double *t_CIA, *lambda_CIA, **opac_CIA;
  int i, j, k, a, b, dum;
  char filename[65];

  FILE *f1;

  /* Allocate Memory */

  opac_CIA_H_el = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_He_H = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_CH4_CH4 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_H2_He = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_H2_CH4 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_H2_H = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_H2_H2 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);
  opac_CIA_CO2_CO2 = dmatrix(0, NTEMP-1, 0, NLAMBDA-1);

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

  /* Fill in CH4 opacities */
  
  opacCH4.T = dvector(0, NTEMP-1);
  opacCH4.P = dvector(0, NPRESSURE-1);
  opacCH4.Plog10 = dvector(0, NPRESSURE-1);
  opacCH4.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacCH4.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacCH4.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacCH4.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacCH4.abundance[j][k] = chem.CH4[j][k];
    }
  }

  strcpy(filename, CH4_FILE);

  ReadOpacTable(opacCH4, filename);
 
  printf("Read CH4 Opacity done\n");

  /* Fill in CO opacities */
  
  opacCO.T = dvector(0, NTEMP-1);
  opacCO.P = dvector(0, NPRESSURE-1);
  opacCO.Plog10 = dvector(0, NPRESSURE-1);
  opacCO.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacCO.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacCO.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacCO.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);

  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacCO.abundance[j][k] = chem.CO[j][k];
    }
  }

  strcpy(filename, CO_FILE);

  ReadOpacTable(opacCO, filename);
 
  printf("Read CO Opacity done\n");

  /* Fill in CO2 opacities */
  
  opacCO2.T = dvector(0, NTEMP-1);
  opacCO2.P = dvector(0, NPRESSURE-1);
  opacCO2.Plog10 = dvector(0, NPRESSURE-1);
  opacCO2.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacCO2.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacCO2.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacCO2.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);

  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){ 
      opacCO2.abundance[j][k] = chem.CO2[j][k];
    }
  }

  strcpy(filename, CO2_FILE);

  ReadOpacTable(opacCO2, filename);
 
  printf("Read CO2 Opacity done\n");

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

  /* Fill in NH3 opacities */
  
  opacNH3.T = dvector(0, NTEMP-1);
  opacNH3.P = dvector(0, NPRESSURE-1);
  opacNH3.Plog10 = dvector(0, NPRESSURE-1);
  opacNH3.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacNH3.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacNH3.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacNH3.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);  

  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacNH3.abundance[j][k] = chem.NH3[j][k];
    }
  }

  strcpy(filename, NH3_FILE);

  ReadOpacTable(opacNH3, filename);
 
  printf("Read NH3 Opacity done\n");

  /* Fill in O2 opacities */
  
  opacO2.T = dvector(0, NTEMP-1);
  opacO2.P = dvector(0, NPRESSURE-1);
  opacO2.Plog10 = dvector(0, NPRESSURE-1);
  opacO2.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacO2.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacO2.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacO2.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);

  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacO2.abundance[j][k] = chem.O2[j][k];
    }
  }
  
  strcpy(filename, O2_FILE);

  ReadOpacTable(opacO2, filename);
 
  printf("Read O2 Opacity done\n");

  /* Fill in O3 opacities */
  
  opacO3.T = dvector(0, NTEMP-1);
  opacO3.P = dvector(0, NPRESSURE-1);
  opacO3.Plog10 = dvector(0, NPRESSURE-1);
  opacO3.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacO3.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacO3.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacO3.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);

  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacO3.abundance[j][k] = chem.O3[j][k];
    }
  }

  strcpy(filename, O3_FILE);

  ReadOpacTable(opacO3, filename);
 
  printf("Read O3 Opacity done\n");

  /* Fill in scattering coefficients */
  
  opacscat.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacscat.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacscat.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }



  /* Fill in collision-induced opacities */


  f1 = fopen("DATA/SET_3/opacCIA_hires.dat", "r");
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
      fscanf(f1, "%le %le %le %le %le %le %le %le %le", &atmos.lambda[i],
                                                        &opac_CIA_He_H[k][i],
                                                        &opac_CIA_CH4_CH4[k][i],
                                                        &opac_CIA_H2_He[k][i],
                                                        &opac_CIA_H2_CH4[k][i],
                                                        &opac_CIA_H2_H[k][i],
                                                        &opac_CIA_H2_H2[k][i],
                                                        &opac_CIA_CO2_CO2[k][i],
                                                        &opac_CIA_H_el[k][i]);
    }
  }
 
  opacCIA.kappa = calloc(NLAMBDA, sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacCIA.kappa[i] = calloc(NPRESSURE, sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacCIA.kappa[i][j] = calloc(NTEMP, sizeof(double));
    }
  }

  for (i=0; i<NLAMBDA; i++){
    for (j=0; j<NPRESSURE; j++){
      for (k=0; k<NTEMP; k++){

      opacCIA.kappa[i][j][k] += opac_CIA_H_el[k][i] *
        (chem.H[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k])) *
        (chem.el[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]));

      opacCIA.kappa[i][j][k] += opac_CIA_He_H[k][i] *
        chem.He[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]) *
        chem.H[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]);


      opacCIA.kappa[i][j][k] += opac_CIA_CH4_CH4[k][i] *
        chem.CH4[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]) *
        chem.CH4[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]);


      opacCIA.kappa[i][j][k] += opac_CIA_H2_He[k][i] *
        chem.H2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]) *
        chem.He[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]);


      opacCIA.kappa[i][j][k] += opac_CIA_H2_CH4[k][i] *
        chem.H2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]) *
        chem.CH4[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]);


      opacCIA.kappa[i][j][k] += opac_CIA_H2_H[k][i] *
        chem.H2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]) *
        chem.H[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]);


      opacCIA.kappa[i][j][k] += opac_CIA_H2_H2[k][i] *
        chem.H2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]) *
        chem.H2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]);


      opacCIA.kappa[i][j][k] += opac_CIA_CO2_CO2[k][i] *
        chem.CO2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]) *
        chem.CO2[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]);
      }
    }
  }

  fclose(f1);
	
/* Rayleigh scattering */

  /* (Polarizabilities from the CRC Handbook) */

  for (i=0; i<NLAMBDA; i++) {
    for (j=0; j<NPRESSURE; j++) {
      for (k=0; k<NTEMP; k++) {
	  opacscat.kappa[i][j][k] +=
	  (8.0*PI/3.0) * SQ(0.80e-30) *
	  SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	  chem.H2[j][k]*chem.P[j] / (KBOLTZMANN * chem.T[k])
	  +
	  (8.0*PI/3.0) * SQ(0.21e-30) *
	  SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	  chem.He[j][k]*chem.P[j] / (KBOLTZMANN * chem.T[k])
	  +
	  (8.0*PI/3.0) * SQ(1.74e-30) *
	  SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	  chem.N2[j][k]*chem.P[j] / (KBOLTZMANN * chem.T[k])
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
	  chem.CH4[j][k]*chem.P[j] / (KBOLTZMANN * chem.T[k])
	  +
	  (8.0*PI/3.0) * SQ(1.58e-30) *
	  SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	  chem.O2[j][k]*chem.P[j] / (KBOLTZMANN * chem.T[k])
	  +
	  (8.0*PI/3.0) * SQ(3.21e-30) *
	  SQ(2.0*PI/ atmos.lambda[i]) * SQ(2.0*PI/ atmos.lambda[i]) *
	  chem.O3[j][k]*chem.P[j] / (KBOLTZMANN * chem.T[k]);
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
    opac.T[k] = opacCO2.T[k];

  for (j=0; j<NPRESSURE; j++) {
    opac.P[j] = opacCO2.P[j];
    opac.Plog10[j] = opacCO2.Plog10[j];
  }
 
  for (i=0; i<NLAMBDA; i++) {
    for (j=0; j<NPRESSURE; j++) {
      for (k=0; k<NTEMP; k++) {
          opac.kappa[i][j][k] = opacCH4.kappa[i][j][k] + opacCO.kappa[i][j][k]
	                    + opacCO2.kappa[i][j][k] + opacH2O.kappa[i][j][k]
	                    + opacNH3.kappa[i][j][k] + opacO2.kappa[i][j][k]
 	                    + opacO3.kappa[i][j][k] + opacscat.kappa[i][j][k]
 	                    + opacCIA.kappa[i][j][k];
      }
    }
  }

  
  /* Free uneeded opacity structures and chemistry table */

  FreeOpacTable(opacCH4);
  FreeOpacTable(opacCO);
  FreeOpacTable(opacCO2);
  FreeOpacTable(opacH2O);
  FreeOpacTable(opacNH3);
  FreeOpacTable(opacO2);
  FreeOpacTable(opacO3);
  FreeOpacTable(opacCIA);
  FreeOpacTable(opacscat);
  FreeChemTable();

  free_dmatrix(opac_CIA_H_el, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_He_H, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_CH4_CH4, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_H2_He, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_H2_CH4, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_H2_H, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_H2_H2, 0, NTEMP-1, 0, NLAMBDA-1);
  free_dmatrix(opac_CIA_CO2_CO2, 0, NTEMP-1, 0, NLAMBDA-1);
}

/* ------- end -------------- TotalOpac.c ------------------------ */
