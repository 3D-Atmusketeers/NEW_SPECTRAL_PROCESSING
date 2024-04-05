/*--------------------------- lin_interpol.c ---------------------

Author: Eliza Miller-Ricci (emillerricci@cfa.harvard.edu)
Last modified:

------------------------------------------------------------------ */

#include <stdio.h>
#include <math.h>
#include "rpctypes.h"

/* --- Global variables ------------------------------------------ */

/* --- Function prototypes --------------------------------------- */

void Locate(int n, double *array, double value, int *ilow);
double lint(double x1, double y1, double x2, double y2, double x);
double lint2D(double x1, double x2, double y1, double y2, double z1, 
	      double z2, double z3, double z4, double x, double y);

/* ------- begin ---------------- lint3D.c ----------------------- */

/* 3-D Linear interpolation procedure
 
------------------------------------------------------------------ */

double lint3D(double x1, double x2, double y1, double y2, double z1, 
	      double z2, double f1, double f2, double f3, double f4, 
	      double f5, double f6, double f7, double f8, double x, 
	      double y, double z)
{
  double f_int1, f_int2, f;
  
  f_int1 = lint2D(x1, x2, y1, y2, f1, f2, f3, f4, x, y);
  f_int2 = lint2D(x1, x2, y1, y2, f5, f6, f7, f8, x, y);
  f = lint(z1, f_int1, z2, f_int2, z);

  return f;
}

/* ------- end ------------------ lint3D.c ----------------------- */

/* ------- begin ---------------- lint2D.c ----------------------- */

/* 1-D Linear interpolation procedure:  Finds the point on a 
   plane (x,y,z) that lies in this square: 
   * (x1,y1,z1)   * (x2,y1,z2)

   * (x1,y2,z3)   * (x2,y2,z4)
------------------------------------------------------------------ */

/* double lint2D(double xx[], double yy[], double zz[][], int nx, int ny,  */
/* 	      double x, double y) */
/* {  */
/*   double z; */
/*   int a, b; */
  
/*   locate(xx, nx-1, x, &a);   */
/*   z = 0.0; */
/*   return z; */
/* } */
double lint2D(double x1, double x2, double y1, double y2, double z1, 
	      double z2, double z3, double z4, double x, double y)
{
  double z_int1, z_int2, z;
  
  z_int1 = lint(x1, z1, x2, z2, x);
  z_int2 = lint(x1, z3, x2, z4, x);
  z = lint(y1, z_int1, y2, z_int2, y);

  return z;
}

/* ------- end ------------------ lint2D.c ----------------------- */

/* ------- begin ---------------- lint.c ------------------------- */

/* 1-D Linear interpolation procedure:  Finds the point on a 
   line (x,y) that lies between (x1,y1) and (x2,y2).
------------------------------------------------------------------ */

double lint(double x1, double y1, double x2, double y2, double x)
{
  double slope, intcept, y;
  slope = (y2 - y1) / (x2 - x1);
  intcept = y1 - slope * x1;
  y = intcept + slope * x;
  return y;
}

/* ------- end ------------------ lint.c ------------------------- */

/* ------- begin ---------------- Linear.c ----------------------- */

double Linear(int n, double xx[], double yy[], double x)
{
  bool_t ascend;
  double xmin, xmax, slope, intcept, y;
  int i; 
  
  ascend = (xx[1] > xx[0]) ? TRUE : FALSE;
  xmin = (ascend) ? xx[0] : xx[n-1];
  xmax = (ascend) ? xx[n-1] : xx[0];

  if (x <= xmin){
    y = (ascend) ? yy[0] : yy[n-1];
    /* printf("WARNING: Linear: Interpolated value off the grid\n"); */
  }  
  else if (x >= xmax){
    y = (ascend) ? yy[n-1] : yy[0];
    /* printf("WARNING: Linear: Interpolated value off the grid\n"); */
  }
  else {
    Locate(n, xx, x, &i);
    slope = (yy[i+1] - yy[i]) / (xx[i+1] - xx[i]);
    intcept = yy[i] - slope * xx[i];
    y = intcept + slope * x;
  }
  return y;
  

}

/* ------- end ------------------ Linear.c ----------------------- */

/* ------- begin ---------------- Map1.c ------------------------- */

void Map1(double *xold, double *fold, int nold, double *xnew, double *fnew, int nnew)
{
  int l, ll, k, lm1, lm2, lp1;
  double ab, af, bb, bf, cb, cf, d, db, df, wt;

  l = 2;
  ll = 0;

  for(k=1; k<=nnew; k++){

    do{
      l = l + 1;
    } while(l <= nold && xnew[k] >= xold[l]);

    if (l > nold) 
      l = nold;

    if (l > 2 && l < nold){

/* PARABOLIC CASE */

      if (l != ll){

	if (l > 3 && l == ll+1){
	  ab = af;
	  bb = bf;
	  cb = cf;
	  
	}else{

/* MUST COMPUTE THE BACKWARD COEFFICIENTS */

	  lm1 = l - 1;
	  lm2 = l - 2;
	  d = (fold[lm1-1] - fold[lm2-1]) / (xold[lm1-1] - xold[lm2-1]);
	  cb = ((fold[l-1] - fold[lm1-1]) / (xold[l-1] - xold[lm1-1]) - d) / 
	    (xold[l-1] - xold[lm2-1]);
	  bb = d + cb * (xold[lm1-1] - xold[lm2-1]);
	  ab = fold[lm1-1];
	}

	lp1 = l + 1;
	lm1 = l - 1;
	d = (fold[l-1] - fold[lm1-1]) / (xold[l-1] - xold[lm1-1]);
	cf = ((fold[lp1-1] - fold[l-1]) / (xold[lp1-1] - xold[l-1]) - d) / 
	  (xold[lp1-1] - xold[lm1-1]);
	bf = d + cf * (xold[l-1] - xold[lm1-1]);
	af = fold[l-1];
	wt = 0.0;
	if (cf != 0.0) 
	  wt = fabs(cf) / (fabs(cf) + fabs(cb));
	ll = l;
      }

      df = xnew[k-1] - xold[l-1];
      db = xnew[k-1] - xold[lm1-1];
      fnew[k-1] = (1.0 - wt) * (af + (bf + cf * df) * df) + wt * 
	(ab + (bb + cb * db) * db);

    }else{

      if (l != ll){
	ll = l;
	lm1 = l - 1;
	af = fold[lm1-1];
	bf = (fold[l-1] - fold[lm1-1]) / (xold[l-1] - xold[lm1-1]);
      }

      fnew[k-1] = af + bf * (xnew[k-1] - xold[lm1-1]);
    }

  }
	
  return;
	
}


/*
double smooth_value(double wav_qext[500][500], int pressure_index, double wavelength_microns, double wavelength_array[500], int max_index) {
    // Ensure max_index does not exceed array bounds
    if (max_index > 500) max_index = 500;

    int index;
    Locate(max_index, wavelength_array, wavelength_microns, &index);

    if (index >= max_index - 1) {
        return wav_qext[pressure_index][index];
    } else {
        // Find the exact position of wavelength_microns between the two closest points
        double lower_wavelength = wavelength_array[index];
        double upper_wavelength = wavelength_array[index + 1];
        double weight = (upper_wavelength - wavelength_microns) / (upper_wavelength - lower_wavelength);
        double next_weight = 1.0 - weight;
        return wav_qext[pressure_index][index] * weight + 
               wav_qext[pressure_index][index + 1] * next_weight;
    }
}
*/

/* ------- end ------------------ Map1.c ------------------------- */
