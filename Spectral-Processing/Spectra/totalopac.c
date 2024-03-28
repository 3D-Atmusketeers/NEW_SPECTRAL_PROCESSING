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
struct Opac opacC2H2;
struct Opac opacAlO;
struct Opac opacNa;
struct Opac opacCa;
struct Opac opacO3;
struct Opac opacK;
struct Opac opacFe;
struct Opac opacCaH;
struct Opac opacAlH;
struct Opac opacH2O;
struct Opac opacVO;
struct Opac opacCO2;
struct Opac opacSH;
struct Opac opacMgH;
struct Opac opacHCN;
struct Opac opacSiO;
struct Opac opacCH4;
struct Opac opacNO;
struct Opac opacNH3;
struct Opac opacH2S;
struct Opac opacHF;
struct Opac opacCO;
struct Opac opacSO2;
struct Opac opacSiH;
struct Opac opacCrH;
struct Opac opacOH;
struct Opac opacC;
struct Opac opacCaO;
struct Opac opacTiO;
struct Opac opacPH3;
struct Opac opacO2;
struct Opac opacNaH;
struct Opac opacFe_plus;
struct Opac opacSiH4;
struct Opac opacMg;
struct Opac opacCa_plus;

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
  //printf("ReadChemTable done\n");

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
  /* Fill in AlO opacities */
  opacAlO.T = dvector(0, NTEMP-1);
  opacAlO.P = dvector(0, NPRESSURE-1);
  opacAlO.Plog10 = dvector(0, NPRESSURE-1);
  opacAlO.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacAlO.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacAlO.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacAlO.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacAlO.abundance[j][k] = chem.AlO[j][k];
    }
  }
  strcpy(filename, AlO_FILE);
  ReadOpacTable(opacAlO, filename);
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
  /* Fill in Ca opacities */
  opacCa.T = dvector(0, NTEMP-1);
  opacCa.P = dvector(0, NPRESSURE-1);
  opacCa.Plog10 = dvector(0, NPRESSURE-1);
  opacCa.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacCa.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacCa.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacCa.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacCa.abundance[j][k] = chem.Ca[j][k];
    }
  }
  strcpy(filename, Ca_FILE);
  ReadOpacTable(opacCa, filename);
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
  /* Fill in Fe opacities */
  opacFe.T = dvector(0, NTEMP-1);
  opacFe.P = dvector(0, NPRESSURE-1);
  opacFe.Plog10 = dvector(0, NPRESSURE-1);
  opacFe.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacFe.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacFe.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacFe.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacFe.abundance[j][k] = chem.Fe[j][k];
    }
  }
  strcpy(filename, Fe_FILE);
  ReadOpacTable(opacFe, filename);
  /* Fill in CaH opacities */
  opacCaH.T = dvector(0, NTEMP-1);
  opacCaH.P = dvector(0, NPRESSURE-1);
  opacCaH.Plog10 = dvector(0, NPRESSURE-1);
  opacCaH.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacCaH.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacCaH.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacCaH.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacCaH.abundance[j][k] = chem.CaH[j][k];
    }
  }
  strcpy(filename, CaH_FILE);
  ReadOpacTable(opacCaH, filename);
  /* Fill in AlH opacities */
  opacAlH.T = dvector(0, NTEMP-1);
  opacAlH.P = dvector(0, NPRESSURE-1);
  opacAlH.Plog10 = dvector(0, NPRESSURE-1);
  opacAlH.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacAlH.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacAlH.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacAlH.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacAlH.abundance[j][k] = chem.AlH[j][k];
    }
  }
  strcpy(filename, AlH_FILE);
  ReadOpacTable(opacAlH, filename);
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
  /* Fill in SH opacities */
  opacSH.T = dvector(0, NTEMP-1);
  opacSH.P = dvector(0, NPRESSURE-1);
  opacSH.Plog10 = dvector(0, NPRESSURE-1);
  opacSH.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacSH.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacSH.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacSH.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacSH.abundance[j][k] = chem.SH[j][k];
    }
  }
  strcpy(filename, SH_FILE);
  ReadOpacTable(opacSH, filename);
  /* Fill in MgH opacities */
  opacMgH.T = dvector(0, NTEMP-1);
  opacMgH.P = dvector(0, NPRESSURE-1);
  opacMgH.Plog10 = dvector(0, NPRESSURE-1);
  opacMgH.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacMgH.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacMgH.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacMgH.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacMgH.abundance[j][k] = chem.MgH[j][k];
    }
  }
  strcpy(filename, MgH_FILE);
  ReadOpacTable(opacMgH, filename);
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
  /* Fill in SiO opacities */
  opacSiO.T = dvector(0, NTEMP-1);
  opacSiO.P = dvector(0, NPRESSURE-1);
  opacSiO.Plog10 = dvector(0, NPRESSURE-1);
  opacSiO.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacSiO.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacSiO.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacSiO.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacSiO.abundance[j][k] = chem.SiO[j][k];
    }
  }
  strcpy(filename, SiO_FILE);
  ReadOpacTable(opacSiO, filename);
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
  /* Fill in NO opacities */
  opacNO.T = dvector(0, NTEMP-1);
  opacNO.P = dvector(0, NPRESSURE-1);
  opacNO.Plog10 = dvector(0, NPRESSURE-1);
  opacNO.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacNO.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacNO.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacNO.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacNO.abundance[j][k] = chem.NO[j][k];
    }
  }
  strcpy(filename, NO_FILE);
  ReadOpacTable(opacNO, filename);
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
  /* Fill in HF opacities */
  opacHF.T = dvector(0, NTEMP-1);
  opacHF.P = dvector(0, NPRESSURE-1);
  opacHF.Plog10 = dvector(0, NPRESSURE-1);
  opacHF.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacHF.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacHF.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacHF.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacHF.abundance[j][k] = chem.HF[j][k];
    }
  }
  strcpy(filename, HF_FILE);
  ReadOpacTable(opacHF, filename);
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
  /* Fill in SO2 opacities */
  opacSO2.T = dvector(0, NTEMP-1);
  opacSO2.P = dvector(0, NPRESSURE-1);
  opacSO2.Plog10 = dvector(0, NPRESSURE-1);
  opacSO2.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacSO2.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacSO2.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacSO2.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacSO2.abundance[j][k] = chem.SO2[j][k];
    }
  }
  strcpy(filename, SO2_FILE);
  ReadOpacTable(opacSO2, filename);
  /* Fill in SiH opacities */
  opacSiH.T = dvector(0, NTEMP-1);
  opacSiH.P = dvector(0, NPRESSURE-1);
  opacSiH.Plog10 = dvector(0, NPRESSURE-1);
  opacSiH.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacSiH.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacSiH.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacSiH.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacSiH.abundance[j][k] = chem.SiH[j][k];
    }
  }
  strcpy(filename, SiH_FILE);
  ReadOpacTable(opacSiH, filename);
  /* Fill in CrH opacities */
  opacCrH.T = dvector(0, NTEMP-1);
  opacCrH.P = dvector(0, NPRESSURE-1);
  opacCrH.Plog10 = dvector(0, NPRESSURE-1);
  opacCrH.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacCrH.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacCrH.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacCrH.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacCrH.abundance[j][k] = chem.CrH[j][k];
    }
  }
  strcpy(filename, CrH_FILE);
  ReadOpacTable(opacCrH, filename);
  /* Fill in OH opacities */
  opacOH.T = dvector(0, NTEMP-1);
  opacOH.P = dvector(0, NPRESSURE-1);
  opacOH.Plog10 = dvector(0, NPRESSURE-1);
  opacOH.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacOH.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacOH.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacOH.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacOH.abundance[j][k] = chem.OH[j][k];
    }
  }
  strcpy(filename, OH_FILE);
  ReadOpacTable(opacOH, filename);
  /* Fill in C opacities */
  opacC.T = dvector(0, NTEMP-1);
  opacC.P = dvector(0, NPRESSURE-1);
  opacC.Plog10 = dvector(0, NPRESSURE-1);
  opacC.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacC.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacC.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacC.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacC.abundance[j][k] = chem.C[j][k];
    }
  }
  strcpy(filename, C_FILE);
  ReadOpacTable(opacC, filename);
  /* Fill in CaO opacities */
  opacCaO.T = dvector(0, NTEMP-1);
  opacCaO.P = dvector(0, NPRESSURE-1);
  opacCaO.Plog10 = dvector(0, NPRESSURE-1);
  opacCaO.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacCaO.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacCaO.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacCaO.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacCaO.abundance[j][k] = chem.CaO[j][k];
    }
  }
  strcpy(filename, CaO_FILE);
  ReadOpacTable(opacCaO, filename);
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
  /* Fill in PH3 opacities */
  opacPH3.T = dvector(0, NTEMP-1);
  opacPH3.P = dvector(0, NPRESSURE-1);
  opacPH3.Plog10 = dvector(0, NPRESSURE-1);
  opacPH3.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacPH3.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacPH3.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacPH3.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacPH3.abundance[j][k] = chem.PH3[j][k];
    }
  }
  strcpy(filename, PH3_FILE);
  ReadOpacTable(opacPH3, filename);
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
  /* Fill in NaH opacities */
  opacNaH.T = dvector(0, NTEMP-1);
  opacNaH.P = dvector(0, NPRESSURE-1);
  opacNaH.Plog10 = dvector(0, NPRESSURE-1);
  opacNaH.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacNaH.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacNaH.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacNaH.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacNaH.abundance[j][k] = chem.NaH[j][k];
    }
  }
  strcpy(filename, NaH_FILE);
  ReadOpacTable(opacNaH, filename);
  /* Fill in Fe_plus opacities */
  opacFe_plus.T = dvector(0, NTEMP-1);
  opacFe_plus.P = dvector(0, NPRESSURE-1);
  opacFe_plus.Plog10 = dvector(0, NPRESSURE-1);
  opacFe_plus.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacFe_plus.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacFe_plus.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacFe_plus.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacFe_plus.abundance[j][k] = chem.Fe_plus[j][k];
    }
  }
  strcpy(filename, Fe_plus_FILE);
  ReadOpacTable(opacFe_plus, filename);
  /* Fill in SiH4 opacities */
  opacSiH4.T = dvector(0, NTEMP-1);
  opacSiH4.P = dvector(0, NPRESSURE-1);
  opacSiH4.Plog10 = dvector(0, NPRESSURE-1);
  opacSiH4.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacSiH4.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacSiH4.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacSiH4.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacSiH4.abundance[j][k] = chem.SiH4[j][k];
    }
  }
  strcpy(filename, SiH4_FILE);
  ReadOpacTable(opacSiH4, filename);
  /* Fill in Mg opacities */
  opacMg.T = dvector(0, NTEMP-1);
  opacMg.P = dvector(0, NPRESSURE-1);
  opacMg.Plog10 = dvector(0, NPRESSURE-1);
  opacMg.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacMg.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacMg.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacMg.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacMg.abundance[j][k] = chem.Mg[j][k];
    }
  }
  strcpy(filename, Mg_FILE);
  ReadOpacTable(opacMg, filename);
  /* Fill in Ca_plus opacities */
  opacCa_plus.T = dvector(0, NTEMP-1);
  opacCa_plus.P = dvector(0, NPRESSURE-1);
  opacCa_plus.Plog10 = dvector(0, NPRESSURE-1);
  opacCa_plus.kappa = malloc(NLAMBDA*sizeof(double));
  for(i=0; i<NLAMBDA; i++){
    opacCa_plus.kappa[i] = malloc(NPRESSURE*sizeof(double));
    for(j=0; j<NPRESSURE; j++){
      opacCa_plus.kappa[i][j] = malloc(NTEMP*sizeof(double));
    }
  }
  opacCa_plus.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  for(j=0; j<NPRESSURE; j++){
    for(k=0; k<NTEMP; k++){
      opacCa_plus.abundance[j][k] = chem.Ca_plus[j][k];
    }
  }
  strcpy(filename, Ca_plus_FILE);
  ReadOpacTable(opacCa_plus, filename);

  printf("=====================================================\n");
  printf("\n==== Debugging H2O Opacities ====\n");
  printf("\tAt lambda index 0, pressure index 0, and temperature index 0:\n");
  printf("\tTemperature: %lf K\n", opacH2O.T[0]);
  printf("\tPressure: %le Pa\n", opacH2O.P[0]);
  printf("\tAbundance: %le\n", opacH2O.abundance[0][0]);
  printf("=====================================================n\n");


  /* Fill in collision-induced opacities */
  //f1 = fopen(CIA_FILE, "r");
  f1 = fopen("DATA/Low-Res/opacCIA.dat", "r");

  if(f1 == NULL){
    printf("\n totalopac.c:\nError opening CIA file: -- No such file or directory\n\n");
  }
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


  printf("\n==== Debug Information for Collisional Opacities ====\n");
  printf("Lambda Range: %e to %e\n", atmos.lambda[0], atmos.lambda[NLAMBDA-1]);
  printf("First 4 Lambdas: %e, %e, %e, %e\n", atmos.lambda[0], atmos.lambda[1], atmos.lambda[2], atmos.lambda[3]);
  printf("Hel Opacity for First Temperature at First 4 Lambdas:\n");
  printf("\tLambda 1: %0.7e\n", opac_CIA_Hel[0][0]);
  printf("\tLambda 2: %0.7e\n", opac_CIA_Hel[0][1]);
  printf("\tLambda 3: %0.7e\n", opac_CIA_Hel[0][2]);
  printf("\tLambda 4: %0.7e\n", opac_CIA_Hel[0][3]);
  printf("\tLambda 5: %0.7e\n", opac_CIA_Hel[0][4]);
  printf("=====================================================\n\n");

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
          opac.kappa[i][j][k] =                              + opacCa_plus.kappa[i][j][k]
                              + opacMg.kappa[i][j][k]
                              + opacSiH4.kappa[i][j][k]
                              + opacFe_plus.kappa[i][j][k]
                              + opacNaH.kappa[i][j][k]
                              + opacO2.kappa[i][j][k]
                              + opacPH3.kappa[i][j][k]
                              + opacTiO.kappa[i][j][k]
                              + opacCaO.kappa[i][j][k]
                              + opacC.kappa[i][j][k]
                              + opacOH.kappa[i][j][k]
                              + opacCrH.kappa[i][j][k]
                              + opacSiH.kappa[i][j][k]
                              + opacSO2.kappa[i][j][k]
                              + opacCO.kappa[i][j][k]
                              + opacHF.kappa[i][j][k]
                              + opacH2S.kappa[i][j][k]
                              + opacNH3.kappa[i][j][k]
                              + opacNO.kappa[i][j][k]
                              + opacCH4.kappa[i][j][k]
                              + opacSiO.kappa[i][j][k]
                              + opacHCN.kappa[i][j][k]
                              + opacMgH.kappa[i][j][k]
                              + opacSH.kappa[i][j][k]
                              + opacCO2.kappa[i][j][k]
                              + opacVO.kappa[i][j][k]
                              + opacH2O.kappa[i][j][k]
                              + opacAlH.kappa[i][j][k]
                              + opacCaH.kappa[i][j][k]
                              + opacFe.kappa[i][j][k]
                              + opacK.kappa[i][j][k]
                              + opacO3.kappa[i][j][k]
                              + opacCa.kappa[i][j][k]
                              + opacNa.kappa[i][j][k]
                              + opacAlO.kappa[i][j][k]
                              + opacC2H2.kappa[i][j][k]
                              + opacscat.kappa[i][j][k]
                              + opacCIA.kappa[i][j][k];
      }
    }
  }


  
  // Free not needed opacity structures and chemistry table
  FreeOpacTable(opacC2H2);
  FreeOpacTable(opacAlO);
  FreeOpacTable(opacNa);
  FreeOpacTable(opacCa);
  FreeOpacTable(opacO3);
  FreeOpacTable(opacK);
  FreeOpacTable(opacFe);
  FreeOpacTable(opacCaH);
  FreeOpacTable(opacAlH);
  FreeOpacTable(opacH2O);
  FreeOpacTable(opacVO);
  FreeOpacTable(opacCO2);
  FreeOpacTable(opacSH);
  FreeOpacTable(opacMgH);
  FreeOpacTable(opacHCN);
  FreeOpacTable(opacSiO);
  FreeOpacTable(opacCH4);
  FreeOpacTable(opacNO);
  FreeOpacTable(opacNH3);
  FreeOpacTable(opacH2S);
  FreeOpacTable(opacHF);
  FreeOpacTable(opacCO);
  FreeOpacTable(opacSO2);
  FreeOpacTable(opacSiH);
  FreeOpacTable(opacCrH);
  FreeOpacTable(opacOH);
  FreeOpacTable(opacC);
  FreeOpacTable(opacCaO);
  FreeOpacTable(opacTiO);
  FreeOpacTable(opacPH3);
  FreeOpacTable(opacO2);
  FreeOpacTable(opacNaH);
  FreeOpacTable(opacFe_plus);
  FreeOpacTable(opacSiH4);
  FreeOpacTable(opacMg);
  FreeOpacTable(opacCa_plus);


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
