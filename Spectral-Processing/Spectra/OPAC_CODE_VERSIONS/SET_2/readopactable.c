
/*----------------------- readopactable.c ------------------------

Author: Sara Seager
Modified by: Eliza Miller-Ricci (emillerricci@cfa.harvard.edu)
Last modified: October 20, 2009

------------------------------------------------------------------ */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "input.h"
#include "atmos.h"
#include "opac.h"
#include "constant.h" 

extern struct Atmos atmos;
 
/* ---------------------------------------------------------------
 * Read in opacity files: kappa(pressure, temperature, lambda)
 * Opacity file actually contains cross section in units of m^2.  
 * Need to multiply by the number density to get kappa in units of
 * m^-1
 * --------------------------------------------------------------- */

/* ------- begin ------------ ReadOpacTable.c -------------------- */

void ReadOpacTable(struct Opac opac, char *filename) {

  int i, j, k, fmt;
  double junk;

  FILE *f1;

  fmt = FORMAT;

  atmos.lambda = malloc(NLAMBDA*sizeof(double));


  switch(fmt){
  case (1): {
    f1 = fopen(filename,"r");
    if(f1 == NULL){
      printf("\nreadopactable.c:\nError opening file: No such file or directory\n\n");
      exit(1);
    }

    fscanf(f1,"%d %d %d", &opac.NP, &opac.NT, &atmos.Nlambda);
    
    for (k=0; k<opac.NT; k++) {
      fscanf(f1,"%le", &opac.T[k]);
    }
    
    for (j=0; j<opac.NP; j++) {
      fscanf(f1,"%le", &opac.P[j]);
      opac.Plog10[j] = log10(opac.P[j]);
    }
    
    for (i=0; i<NLAMBDA; i++) {
      fscanf(f1,"%le", &atmos.lambda[i]);
    }
    
    for (j=0; j<opac.NP; j++) {
      for (k=0; k<opac.NT; k++) {
	for (i=0; i<NLAMBDA; i++) {
	  fscanf(f1,"%le", &opac.kappa[i][j][k]);
	}
      }
    }
    fclose(f1);
    break;
  }

  case(2): {
    opac.NP = NPRESSURE;
    opac.NT = NTEMP; 
    atmos.Nlambda = NLAMBDA;

    f1 = fopen(filename,"r");
    if(f1 == NULL){
      printf("\nreadopactable.c:\nError opening file: No such file or directory\n\n");
      exit(1);
    }

    for (k=0; k<opac.NT; k++) {
      fscanf(f1,"%le", &opac.T[k]);
    }
    
    for (j=0; j<opac.NP; j++) {
      fscanf(f1,"%le", &opac.P[j]);
      opac.Plog10[j] = log10(opac.P[j]);
    }

    for (i=0; i<NLAMBDA; i++) {
      fscanf(f1,"%le", &atmos.lambda[i]);
      for (j=0; j<opac.NP; j++) {
	fscanf(f1,"%le", &junk);
	for (k=0; k<opac.NT; k++) {
	  fscanf(f1,"%le", &opac.kappa[i][j][k]);
	      
	  /*   kappa in file is actually cross section, sigma.  
	       Need to multiply by number density */
	  
	  opac.kappa[i][j][k] *= opac.abundance[j][k] * opac.P[j] /
	                         (KBOLTZMANN * opac.T[k]);
	}
      }
    }
    
    fclose(f1);
    printf("opac %e %e %e\n", atmos.lambda[NLAMBDA-1], opac.P[0], 
	   opac.T[0]);
    break;
  }

  default: {
    printf("Invalid format for opacity table\n\n");
    exit(1);
  }
    
  }
}

/* ------- end -------------- ReadOpacTable.c -------------------- */

/* ------- begin ------------ FreeOpacTable.c -------------------- */

void FreeOpacTable(struct Opac opac)
{

  free(opac.T);
  free(opac.P);
  free(opac.Plog10);
  free(opac.kappa);
  free(opac.kappa_pl);

}

/* ------- end -------------- FreeOpacTable.c -------------------- */
