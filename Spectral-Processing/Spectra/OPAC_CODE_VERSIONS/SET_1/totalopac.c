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

struct Opac opacC2H2;
struct Opac opacCH4;
struct Opac opacCO;
struct Opac opacCO2;
struct Opac opacFeH;
struct Opac opacH2O;
struct Opac opacH2S;
struct Opac opacHCN;
struct Opac opacK;
struct Opac opacNa;
struct Opac opacNH3;
struct Opac opacTiO;
struct Opac opacVO;
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
         **opac_CIA_N2N2, **opac_CIA_O2CO2, **opac_CIA_O2N2, **opac_CIA_O2O2;
  double *t_CIA, *lambda_CIA, **opac_CIA;
  int i, j, k, a, b, dum;
  char filename[65];

  double m_H2 = 2.0158;
  double m_H = 1.0079;
  double m_He = 4.002602;
  double m_H2O = 18.0152;
  double m_CH4 = 16.0423;
  double m_CO = 28.010;
  double m_CO2 = 44.010;
  double m_O = 15.9994;
  double m_C = 12.0107;
  double m_N = 14.0067;
  double m_NH3 = 17.031;
  double m_N2 = 28.0134;
  double m_Na = 22.988977;
  double m_Naxx  = 39.997;
  double m_K = 39.0983;
  double m_Kxx = 56.1056;
  double m_O2 = 31.9988;
  double m_O3 = 47.9982;

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

  /* Fill in C2H2 opacities */

  opacC2H2.T = dvector(0, NTEMP-1);
  opacC2H2.P = dvector(0, NPRESSURE-1);
  opacC2H2.Plog10 = dvector(0, NPRESSURE-1);
  opacC2H2.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacC2H2.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacC2H2.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacC2H2.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);

  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacC2H2.abundance[j][k] = chem.C2H2[j][k];
    }
  }

  strcpy(filename, C2H2_FILE);

  ReadOpacTable(opacC2H2, filename);

  printf("Read C2H2 Opacity done\n");

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


  /* Fill in FeH opacities */

  opacFeH.T = dvector(0, NTEMP-1);
  opacFeH.P = dvector(0, NPRESSURE-1);
  opacFeH.Plog10 = dvector(0, NPRESSURE-1);
  opacFeH.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacFeH.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacFeH.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacFeH.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);

  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacFeH.abundance[j][k] = chem.FeH[j][k];
    }
  }

  strcpy(filename, FeH_FILE);

  ReadOpacTable(opacFeH, filename);

  printf("Read FeH Opacity done\n");

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


  /* Fill in H2S opacities */

  opacH2S.T = dvector(0, NTEMP-1);
  opacH2S.P = dvector(0, NPRESSURE-1);
  opacH2S.Plog10 = dvector(0, NPRESSURE-1);
  opacH2S.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacH2S.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacH2S.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacH2S.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);

  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacH2S.abundance[j][k] = chem.H2S[j][k];
    }
  }

  strcpy(filename, H2S_FILE);

  ReadOpacTable(opacH2S, filename);

  printf("Read H2S Opacity done\n");



  /* Fill in HCN opacities */

  opacHCN.T = dvector(0, NTEMP-1);
  opacHCN.P = dvector(0, NPRESSURE-1);
  opacHCN.Plog10 = dvector(0, NPRESSURE-1);
  opacHCN.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacHCN.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacHCN.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacHCN.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);

  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacHCN.abundance[j][k] = chem.HCN[j][k];
    }
  }

  strcpy(filename, HCN_FILE);

  ReadOpacTable(opacHCN, filename);

  printf("Read HCN Opacity done\n");


  /* Fill in K opacities */

  opacK.T = dvector(0, NTEMP-1);
  opacK.P = dvector(0, NPRESSURE-1);
  opacK.Plog10 = dvector(0, NPRESSURE-1);
  opacK.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacK.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacK.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacK.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);

  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacK.abundance[j][k] = chem.K[j][k];
    }
  }

  strcpy(filename, K_FILE);

  ReadOpacTable(opacK, filename);

  printf("Read K Opacity done\n");

  /* Fill in Na opacities */

  opacNa.T = dvector(0, NTEMP-1);
  opacNa.P = dvector(0, NPRESSURE-1);
  opacNa.Plog10 = dvector(0, NPRESSURE-1);
  opacNa.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacNa.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacNa.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacNa.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);

  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacNa.abundance[j][k] = chem.Na[j][k];
    }
  }

  strcpy(filename, Na_FILE);

  ReadOpacTable(opacNa, filename);

  printf("Read Na Opacity done\n");


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

  /* Fill in TiO opacities */

  opacTiO.T = dvector(0, NTEMP-1);
  opacTiO.P = dvector(0, NPRESSURE-1);
  opacTiO.Plog10 = dvector(0, NPRESSURE-1);
  opacTiO.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacTiO.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacTiO.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacTiO.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);

  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacTiO.abundance[j][k] = chem.TiO[j][k];
    }
  }

  strcpy(filename, TiO_FILE);

  ReadOpacTable(opacTiO, filename);

  printf("Read TiO Opacity done\n");

  /* Fill in VO opacities */

  opacVO.T = dvector(0, NTEMP-1);
  opacVO.P = dvector(0, NPRESSURE-1);
  opacVO.Plog10 = dvector(0, NPRESSURE-1);
  opacVO.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacVO.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacVO.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacVO.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);

  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacVO.abundance[j][k] = chem.VO[j][k];
    }
  }

  strcpy(filename, VO_FILE);

  ReadOpacTable(opacVO, filename);

  printf("Read VO Opacity done\n");


  /* Fill in scattering coefficients */
  
  opacscat.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacscat.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacscat.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }






  /* Fill in collision-induced opacities */


  f1 = fopen("DATA/SET_1/opacCIA_lotemp_isaac_v2.dat", "r");
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
      fscanf(f1, "%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le", &atmos.lambda[i],
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
  printf("CIA: %e   %e\n", atmos.lambda[NLAMBDA-1], opac_CIA_H2H2[NTEMP-1][NLAMBDA-1]);
 
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


      opacCIA.kappa[i][j][k] += opac_CIA_CH4Ar[k][i] *
        chem.CH4[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]) *
        chem.Ar[j][k] * opacH2O.P[j] / (KBOLTZMANN * opacH2O.T[k]);


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
    opac.T[k] = opacCO2.T[k];

  for (j=0; j<NPRESSURE; j++) {
    opac.P[j] = opacCO2.P[j];
    opac.Plog10[j] = opacCO2.Plog10[j];
  }
 
  for (i=0; i<NLAMBDA; i++) {
    for (j=0; j<NPRESSURE; j++) {
      for (k=0; k<NTEMP; k++) {
          opac.kappa[i][j][k] = opacC2H2.kappa[i][j][k]
                              + opacCH4.kappa[i][j][k]
                              + opacCO.kappa[i][j][k]
	                            + opacCO2.kappa[i][j][k]
	                            + opacFeH.kappa[i][j][k]
	                            + opacH2O.kappa[i][j][k]
	                            + opacH2S.kappa[i][j][k]
	                            + opacHCN.kappa[i][j][k]
	                            + opacK.kappa[i][j][k]
	                            + opacNa.kappa[i][j][k]
	                            + opacNH3.kappa[i][j][k]
                              + opacTiO.kappa[i][j][k]
                              + opacVO.kappa[i][j][k]
 	                            + opacscat.kappa[i][j][k]
 	                            + opacCIA.kappa[i][j][k];
      }
    }
  }


  
  /* Free uneeded opacity structures and chemistry table */

  FreeOpacTable(opacC2H2);
  FreeOpacTable(opacCH4);
  FreeOpacTable(opacCO);
  FreeOpacTable(opacCO2);
  FreeOpacTable(opacFeH);
  FreeOpacTable(opacH2O);
  FreeOpacTable(opacH2S);
  FreeOpacTable(opacHCN);
  FreeOpacTable(opacK);
  FreeOpacTable(opacNa);
  FreeOpacTable(opacNH3);
  FreeOpacTable(opacTiO);
  FreeOpacTable(opacVO);
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

}

/* ------- end -------------- TotalOpac.c ------------------------ */
