/*-------------------------- geometry.c --------------------------

Author: Eliza Miller-Ricci (emillerricci@cfa.harvard.edu)
Last modified:

------------------------------------------------------------------ */

/* Routines dealing with the transmission geometry for calculating
a transmission spectrum. 

------------------------------------------------------------------ */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "input.h"
#include "opac.h"
#include "atmos.h"
#include "constant.h"
#include "include.h"

/* --- Global variables ------------------------------------------ */

extern struct Atmos atmos;
extern struct Opac opac;

/* --- Function prototypes --------------------------------------- */

double lint(double x1, double y1, double x2, double y2, double x);
double Radius(double R_pl, double ds[]);
void PathLength(double ds[], double l_chord[], double dtheta[]);

/* ------- begin ---------------- Radius ------------------------- */

/* Given the planet radius, this determines the radius
for the planet plus atmosphere: R.

------------------------------------------------------------------ */

double Radius(double R_pl, double ds[])
{
  double R;
  int i;
  
  R = R_pl;

  /* i=0; */
/*   do{ */
/*     if (i == NTAU) { */
/*       printf("WARNING: P never reaches 1 bar\n"); */
/*       break; */
/*     } */
/*     R += ds[i]; */
/*     i++; */
/*   }while(atmos.P[i] < 1.0e+05); */

/*   if (i < NTAU) { */
/*     R += lint(atmos.P[i-1], 0, atmos.P[i], ds[i], 1.0e+05); */
/*   } */

  for(i=0; i<NTAU; i++){
    R += ds[i];
  }

  return R;
}

/* ------- end ------------------ Radius ------------------------- */

/* ------- begin ---------------- Radius_3d ---------------------- */

/* Given the planet radius, this determines the radius
for the planet plus atmosphere: R.

------------------------------------------------------------------ */

double Radius_3d(double R_pl, double ds[], double lat)
{
  double R;
  int i;
  
  R = R_pl;

  for(i=0; i<NTAU; i++){
    R += ds[i];
  }

  return R;
}

/* ------- end ------------------ Radius_3d ---------------------- */

/* ------- begin ---------------- PathLength --------------------- */

/* Determines line-of-sight path lengths for rays passing
through the atmosphere, given the radial structure of the
atmosphere.  Also determines d_theta steps for later use in
integrating to obtain the in-transit flux.

------------------------------------------------------------------ */

void PathLength(double ds[], double l_chord[], double dtheta[])
{
  double theta[NTAU];
  double d, R_atm, R;
  int j;

  R = Radius(R_PLANET, ds);

  R_atm = 0.0;
  for(j=0; j<NTAU; j++){
    R_atm += ds[j];
  }
  printf("R_atm %e\n", R_atm);
  
  d = R;
  printf("d %e\n", d);
  for(j=0; j<NTAU; j++){

    d -= ds[j];
    theta[j] = asin(d/R);
    l_chord[j] = 2.0 * R * cos(theta[j]);

    if(j==0)
      dtheta[j] = PI/2.0 - theta[j];
    else
      dtheta[j] = theta[j-1] - theta[j];
  }

  return;
}

/* ------- end ------------------ PathLength --------------------- */

/* ------- begin ---------------- Angles ------------------------- */

/* Determines theta and dtheta values for the flux integral.

------------------------------------------------------------------ */

void Angles(double ds[], double theta[], double dtheta[])
{
  double h, R;
  int j;

  R = Radius(R_PLANET, ds);
  
  h = R;
  printf("R (planet & atmosphere) %e\n", h);

  for(j=0; j<NTAU; j++){

    h -= ds[j];
    theta[j] = asin(h/R_STAR);

    if(j==0)
      dtheta[j] = asin(R/R_STAR) - theta[j];
    else
      dtheta[j] = theta[j-1] - theta[j];

  }

  return;  
}

/* ------- end ------------------ Angles ------------------------- */

/* ------- begin ---------------- Angles3d----------------------- */

/* Determines theta and dtheta values for the flux integral.

------------------------------------------------------------------ */

void Angles3d(double ds[], double theta[], double dtheta[], double lat)
{
  double h, R;
  int j;

  R = Radius(R_PLANET, ds);

  h = R;
  printf("R (planet & atmosphere) %e\n", h);

  for(j=0; j<NTAU; j++){

    h -= ds[j];
    theta[j] = asin(h/R_STAR);

    if(j==0)
      dtheta[j] = asin(R/R_STAR) - theta[j];
    else
      dtheta[j] = theta[j-1] - theta[j];

  }

  return;  
}

/* ------- end ------------------ Angles_3d ---------------------- */

/* ------- begin ---------------- Tau_LOS ------------------------ */

/* Determines line-of-sight optical depths tau_tr as a function 
of wavelength and height in the atmosphere.

------------------------------------------------------------------ */

void Tau_LOS(double **kappa_nu, double **tau_tr, double ds[])
{
  double R, a, b, dl[NTAU], phi;
  int i, j, k;

  R = Radius(R_PLANET, ds);
  printf("R_atmosphere %e\n", R - R_PLANET);

  for(j=0; j<NTAU; j++)
    dl[j] = 0.0;

  a = R;

  for(j=0; j<NTAU; j++){

    a -= ds[j];
    b = R;

    for(k=0; k<=j; k++){
      dl[k] = 2.0 * pow(SQ(b) - SQ(a), 0.5);
      phi = acos(a/b)*180.0/PI;
      if(j == NTAU-1)
	printf("PHI: %d\t%f\n", k, phi+90.0);
      b -= ds[k];
    }

    for(k=0; k<=j; k++){
      if (k != j)
	dl[k] -= dl[k+1];
    }
    
    for(k=0; k<=j; k++){
      for(i=0; i<NLAMBDA; i++){
	
	/* Calculates tau_tr by summing up kappa * dl 
	   for shells along the line of sight */
	tau_tr[i][j] += kappa_nu[i][k] * dl[k];

      }
    }
  }
  
  return;
}

/* ------- end ------------------ Tau_LOS ------------------------ */
