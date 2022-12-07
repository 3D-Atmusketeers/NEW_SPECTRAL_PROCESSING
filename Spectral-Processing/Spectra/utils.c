#include <stdio.h>
#include <math.h>
#include "rpctypes.h"
#include <stdlib.h>
#include "include.h"


/* -------------------------------------- Hunt.c --------------------

       Version:       rh1.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Feb 16 14:40:59 1999 --

       --------------------------                      ----------RH-- */

/* --- Find index of value in array (cf., Num. Recipes, p 91).
       Use previous value of ilow to hunt up or down list and bracket value.
       --                                              -------------- */

/* ------- begin -------------------------- Hunt.c ------------------ */

void Hunt(int n, double *array, double value, int *ilow)
{
  bool_t ascend;
  int    ihigh, index, increment;

  ascend = (array[n-1] > array[0]) ? TRUE : FALSE;
  if ((*ilow <= 0)  ||  (*ilow > n-1)) {

    /* --- Input guess not useful here, go to bisection --  --------- */

    *ilow = 0;
    ihigh = n;
  } else {

    /* --- Else hunt up or down to bracket value --    -------------- */ 

    increment = 1;
    if (((value >= array[*ilow]) ? TRUE : FALSE) == ascend) {
      ihigh = *ilow + increment;
      if (*ilow == n-1) return;

      /* --- Hunt up --                                -------------- */

      while (((value >= array[ihigh]) ? TRUE : FALSE) == ascend) {
	*ilow = ihigh;
	increment += increment;
	ihigh = *ilow + increment;
        if (ihigh >= n) { ihigh = n;  break; }
      }
    } else {
      ihigh = *ilow;
      if (*ilow == 0) return;

      /* --- Hunt down --                              -------------- */

      while (((value <= array[*ilow]) ? TRUE : FALSE) == ascend) {
	ihigh = *ilow;
	increment += increment;
	*ilow = ihigh - increment;
        if (*ilow <= 0) { *ilow = 0;  break; }
      }
    }
  }
  /* --- Bisection algorithm --                        -------------- */

  if (ascend) {
    while (ihigh - *ilow > 1) {
      index = (ihigh + *ilow) >> 1;
      if (value >= array[index])
	*ilow = index;
      else
	ihigh = index;
    }
  } else {
    while (ihigh - *ilow > 1) {
      index = (ihigh + *ilow) >> 1;
      if (value <= array[index])
	*ilow = index;
      else
	ihigh = index;
    }
  }
}
/* ------- end ---------------------------- Hunt.c ------------------ */

/* ---------------------------------------- Locate.c ---------------- */

/* --- Find index of value in array (cf., Num. Recipes, p 90).

 Note: The Num. Recipes routine does not give the correct index
       for values that are exactly equal to an array value!
       --                                              -------------- */
 
/* ------- begin -------------------------- Locate.c ---------------- */

void Locate(int n, double *array, double value, int *ilow)
{
  bool_t ascend;
  int    ihigh, index;

  ascend = (array[n-1] > array[0]) ? TRUE : FALSE;
  *ilow = 0;  ihigh = n;

  if (ascend) {
    while (ihigh - *ilow > 1) {
      index = (ihigh + *ilow) >> 1;
      if (value >= array[index])
	*ilow = index;
      else
	ihigh = index;
    }
  } else {
    while (ihigh - *ilow > 1) {
      index = (ihigh + *ilow) >> 1;
      if (value <= array[index])
	*ilow = index;
      else
	ihigh = index;
    }
  }
}
/* ------- end ---------------------------- Locate.c ---------------- */

/* ---------------------------------------- Deriv.c ----------------- */

/* --- Derivative routine adapted from atlas (atlas9_utils.f)
       Computes derivative of f w/ respect to x and outputs dfdx
       --                                              -------------- */
 
/* ------- begin -------------------------- Deriv.c ----------------- */

void Deriv(double *x, double *f, double *dfdx, int n)
{
  int i;
  double d, d1, s, scale, tan, tan1;

  dfdx[0] = (f[1] - f[0]) / (x[1] - x[0]);
  dfdx[n-1]  = (f[n-1] - f[n-2]) / (x[n-1] - x[n-2]);

  if(n > 2){

    s = fabs(x[1] - x[0]) / (x[1] - x[0]);

    for(i=1; i < n-1; i++){
      scale = MAX(MAX(fabs(f[i-1]), fabs(f[i])), fabs(f[i+1])) / fabs(x[i]);
      if(scale == 0.0)
	scale = 1.0;
      d1 = (f[i+1] - f[i]) / (x[i+1] - x[i]) / scale;
      d =  (f[i] - f[i-1]) / (x[i]   - x[i-1]) / scale;
      tan1 = d1 / (s * sqrt(1.0 + SQ(d1)) + 1.0);
      tan =  d  / (s * sqrt(1.0 + SQ(d)) + 1.0);
      dfdx[i] = (tan1 + tan) / (1.0 - tan1 * tan) * scale;
   
    }
  }

}
/* ------- end ---------------------------- Deriv.c ----------------- */

/* ---------------------------------------- Integ.c ----------------- */

/* --- Numerical integral routine following the left-hand rule
       Computes integral of f dx and outputs fint
       --                                              -------------- */
 
/* ------- begin -------------------------- Integ.c ----------------- */

void Integ(double *x, double *f, double *fint, int n, double start)
{
  double dx[n];
  int i, j;

  if (n < 1){
    printf("Integ Error: Number of elements in array is less than 1\n");
    exit(1);
  }

  for(i=0; i<n; i++){
    if(i==0)
      dx[i] = x[i+1]-x[i];
    else
      dx[i] = x[i] - x[i-1];
  }

  for(i=0; i<n; i++){
    fint[i] = start;
    for(j=0; j<=i; j++){
      fint[i] += f[j]*dx[j];
    }
  }

}

/* ------- end ---------------------------- Integ.c ----------------- */

/* ------- begin -------------------------- Sign.c ------------------ */

/* --- Determines the sign of a double.  Returns +1 if the number is
       positive and -1 if the number is negative.  ------------------ */

int Sign(double n)
{
  if(n >= 0)
    return 1;
  else
    return -1;
}

/* ------- end ---------------------------- Sign.c ------------------ */

/* ---------------------------------------- Expi.c ------------------ */

/* --- Exponential integral for positive arguments after Cody and
       Thatcher, Math. of Comp., 22, 641 (1968) - adapted from 
       atlas9.utils.f                              ------------------ */

/* ------- begin -------------------------- Expi.c ------------------ */

double Expi(int n, double x)
{
  double exp_i, a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, b4, c0, c1, 
    c2, c3, c4, c5, c6, d1, d2, d3, d4, d5, d6, e0, e1, e2, e3, e4,
    e5, e6, f1, f2, f3, f4, f5, f6, x1, ex, ex1;
  int i;

  a0 = -44178.5471728217;
  a1 =  57721.7247139444;
  a2 =  9938.31388962037;
  a3 =  1842.11088668000;
  a4 =  101.093806161906;
  a5 =  5.03416184097568;
  b0 =  76537.3323337614;
  b1 =  32597.1881290275;
  b2 =  6106.10794245759;
  b3 =  635.419418378382;
  b4 =  37.2298352833327;
  c0 =  4.65627107975096e-07;
  c1 =  0.999979577051595;
  c2 =  9.04161556946329;
  c3 =  24.3784088791317;
  c4 =  23.0192559391333;
  c5 =  6.90522522784444;
  c6 =  0.430967839469389;
  d1 =  10.0411643829054;
  d2 =  32.4264210695138;
  d3 =  41.2807841891424;
  d4 =  20.4494785013794;
  d5 =  3.31909213593302;
  d6 =  0.103400130404874;
  e0 = -0.999999999998447;
  e1 = -26.6271060431811;
  e2 = -241.055827097015;
  e3 = -895.927957772937;
  e4 = -1298.85688746484;
  e5 = -545.374158883133;
  e6 = -5.66575206533869;
  f1 =  28.6271060422192;
  f2 =  292.310039388533;
  f3 =  1332.78537748257;
  f4 =  2777.61949509163;
  f5 =  2404.01713225909;
  f6 =  631.657483280800;
  x1 = -1.0e+20;

  if(x != x1){
    x1 = x;
    ex = exp(-x1);

    if(x1 > 4.0){
      ex1 = (ex + ex * (e0 + (e1 + (e2 + (e3 + (e4 + (e5 + e6 / 
	     x1) / x1) / x1) / x1) / x1) / x1) / 
            (x1 + f1 +(f2 + (f3 + (f4 + (f5 + f6 / x1) / x1) / x1) 
	     / x1) / x1)) / x1;

    } else if(x1 > 1.0){
      ex1 = ex * (c6 + (c5 + (c4 + (c3 + (c2 + (c1 + c0 * x1) * 
		  x1) * x1) * x1) * x1) * x1) / (d6 + (d5 + (d4 +
		 (d3 + (d2 + (d1 + x1) * x1) * x1) * x1) * x1) * x1);

    } else if(x1 > 0.0){
      ex1 = (a0 + (a1 + (a2 + (a3 + (a4 + a5 * x1) * x1) * x1) * 
		   x1) * x1) / (b0 + (b1 + (b2 + (b3 + (b4 + x1) * 
		   x1) * x1) * x1) * x1) - log(x1);

    } else {
      ex1 = 0.0;

    }
  }

  exp_i = ex1;

  if(n > 1){
    for(i=1; i<n; i++)
      exp_i = (ex - x1 * exp_i) / i;
  }
  
  return exp_i;
}

/* ------- end ---------------------------- Expi.c ------------------ */
