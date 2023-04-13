#include <stdio.h>

double Planck(double T, double lambda);

void two_stream(int NLAYER, int kmin, double *w0_array, double *g0_array, \
                 const double *temperature_array, const double *tau_array, \
                 double NU, double NU_BIN, double incident_frac, double *dtau_array, \
                 double intensity_vals[])
{
  double mu_1;        // Param for Quadrature or Hemispheric Mean constant
  double mu_0 = 1.0;  // Incident angle of the solar beam
  double mu = 1.0;    // Never used. The angle of the upwards intensity

  // These are indexing values
  int J, L, KINDEX, Z, M;

  // The number of layers in the atmosphere
  int NEW_NLAYER;

  // These are just constants
  double bolz_constant = 1.380649e-23;
  double h_constant    = 6.62607015e-34;

  // These are boundary conditions
  // They vary for long wave and short wave
  double EMIS;
  double RSFX;
  double SFCS_HEMISPHERIC;
  double SFCS_QUADRATURE;

  // The direct radiation from the star
  double DIRECT_HEMISPHERIC[NLAYER - kmin];
  double DIRECT_QUADRATURE[NLAYER - kmin];

  // The flux at the surface of the planet
  // The very top, not the boundary between bottom and atmosphere
  double FLUX_SURFACE_QUADRATURE;
  double FLUX_SURFACE_HEMISPHERIC;

  // Blackbody intensities
  // All SI units, all in the same units as the returned intensity
  double STELLAR_BB;
  double BB_TOP_OF_ATM;
  double BB_BOTTOM_OF_ATM;

  // Scattering and atmosphere parameters
  double W0[NLAYER - kmin];
  double G0[NLAYER - kmin];

  // Scattering and atmosphere parameters
  double TEMPS[NLAYER - kmin];
  double TAUCS[NLAYER - kmin + 1];
  double TAULS[NLAYER - kmin];

  // How to generate the matrix from the toon paper
  double y1[NLAYER - kmin];
  double y2[NLAYER - kmin];
  double y3[NLAYER - kmin];
  double y4[NLAYER - kmin];

  // More Matrix stuff
  double LAMBDAS[NLAYER - kmin];
  double GAMMA[NLAYER - kmin];
  double temp_e_val[NLAYER - kmin];
  double e1[NLAYER - kmin];
  double e2[NLAYER - kmin];
  double e3[NLAYER - kmin];
  double e4[NLAYER - kmin];

  // Matrix Coefficients
  double A[2 * (NLAYER - kmin)];
  double B[2 * (NLAYER - kmin)];
  double D[2 * (NLAYER - kmin)];
  double E[2 * (NLAYER - kmin)];
  
  // Solution Coefficients
  double AS[2 * (NLAYER - kmin)];
  double DS[2 * (NLAYER - kmin)];
  double X[2 * (NLAYER - kmin)];
  double Y[2 * (NLAYER - kmin)];

  // Temporary values, makes math easier for me
  double temp_gamma_val[NLAYER - kmin];

  // Planck Intensities
  // B1 is the slope of intensities between layers
  double B0[NLAYER - kmin];
  double B1[NLAYER - kmin];
  double TEMP[NLAYER - kmin];
  double DTAUS[NLAYER - kmin];

  // Stuff for solving the blackbody equation
  double Bnu, twohnu3_c2, hc_Tkla;
  double temp_val_1, temp_val_2;

  // Matrix values, related to the intensity at each layer
  // These set differently for the quadrature and hemipsheric solutions
  double CP[NLAYER - kmin];
  double CPB[NLAYER - kmin];
  double CM[NLAYER - kmin];
  double CMB[NLAYER - kmin];

  // The Source Function Variables
  double SOURCE_G[NLAYER - kmin];
  double SOURCE_H[NLAYER - kmin];
  double SOURCE_J[NLAYER - kmin];
  double SOURCE_K[NLAYER - kmin];
  double ALPHA_1[NLAYER - kmin];
  double ALPHA_2[NLAYER - kmin];
  double SIGMA_1[NLAYER - kmin];
  double SIGMA_2[NLAYER - kmin];
  double SOURCE_Y1[NLAYER - kmin];
  double SOURCE_Y2[NLAYER - kmin];
  double source_temp[NLAYER - kmin];

  // Upward and downwards intensities
  double INTENSITY_DOWN[NLAYER - kmin];
  double INTENSITY_UP[NLAYER - kmin];
  double TMI[NLAYER - kmin];

  // Final values for long wave radiation
  double HEMISPHERIC_TWO_STREAM[NLAYER - kmin];
  double HEMISPHERIC_SOURCE_FNC[NLAYER - kmin];

  double TOTAL_FLUX;

  // Final values for short wave radiation
  double QUADRATURE_TWO_STREAM[NLAYER - kmin];
  double QUADRATURE_SOURCE_FNC[NLAYER - kmin];

  // Top layer values short wave
  double QUADRATURE_TWO_STREAM_TOP;
  double QUADRATURE_SOURCE_FNC_TOP;

  // The number of layers
  // Sometimes the inputs are bad and the top
  NEW_NLAYER = NLAYER - kmin;

  // Calculate the flux at the top of the atmosphere from the star
  temp_val_1 = (2.0 * h_constant * (NU * NU * NU)) / (CLIGHT * CLIGHT);
  temp_val_2 = exp(h_constant * NU / (bolz_constant * STELLAR_TEMP)) - 1.0;
  STELLAR_BB = temp_val_1 * (1.0 / temp_val_2);

  // Get the flux at the surface of the atmpsphere
  FLUX_SURFACE_HEMISPHERIC = 0;
  FLUX_SURFACE_QUADRATURE = incident_frac * PI * STELLAR_BB * pow((R_STAR / ORB_SEP), 2.0);
  

  // Assign arrays for the characteristics of the atmosphere based on the input arrays

  TAUCS[0] = 0.0;
  for (J=0; J<NEW_NLAYER; J++)
  {
    W0[J] = w0_array[J+kmin];
    G0[J] = g0_array[J+kmin];
    TEMPS[J] = temperature_array[J+kmin];

    //TEMPS[J] = 1000;
    //W0[J] = 0;//W0_VAL;
    //G0[J] = 0;//G0_VAL;

    DTAUS[J]   = dtau_array[J+kmin];
    TAULS[J]   = dtau_array[J+kmin];
    TAUCS[J+1] = TAUCS[J]+TAULS[J];

    DIRECT_QUADRATURE[J]  = mu_0 * PI * FLUX_SURFACE_QUADRATURE * exp(-1.0 * (TAUCS[J] + TAULS[J]) / mu_0);
    DIRECT_HEMISPHERIC[J] = 0.0;
  }

  // I THINK I MAY NEED TO GET RID OF THIS
  TAUCS[0] = TAULS[0];
  
  // Calculate the intensity at the top of the atmosphere
  temp_val_1 = (2.0 * h_constant * (NU * NU * NU)) / (CLIGHT * CLIGHT);
  temp_val_2 = exp(h_constant * NU / (bolz_constant * TEMPS[0])) - 1.0;
  BB_TOP_OF_ATM = temp_val_1 * (1.0 / temp_val_2);

  // Calculate the intensity at the bottom of the atmosphere
  temp_val_1 = (2.0 * h_constant * (NU * NU * NU)) / (CLIGHT * CLIGHT);
  temp_val_2 = exp(h_constant * NU / (bolz_constant * TEMPS[NEW_NLAYER-1])) - 1.0;
  BB_BOTTOM_OF_ATM = temp_val_1 * (1.0 / temp_val_2);


  //***************************************************************************************
  //***************************************************************************************
  //***************************************************************************************
  //***************************************************************************************
  //***************************************************************************************



  //**************************************************
  //*               Long Wave Solution               *
  //**************************************************

  // Boundary conditions
  mu_1 = 0.5;

  EMIS = 1.0;
  RSFX = 1.0;

  SFCS_HEMISPHERIC = EMIS * PI * BB_BOTTOM_OF_ATM;

  // HERE WE FIND LAYER PROPERTIES FOLLOWING GENERAL SCHEME
  // OF MEADOR AND WEAVOR. THEN WE SET UP LAYER PROPERTIES
  // NEEDED FOR MATRIX.

  for(J=0; J<NEW_NLAYER; J++)
  {
    y1[J]    =  (2.0 - (W0[J] * (1.0 + G0[J])));
    y2[J]    =  (W0[J] * (1.0 - G0[J]));

    LAMBDAS[J]    =  sqrt(fabs(pow(y1[J], 2.0) - pow(y2[J], 2.0)));
    GAMMA[J]  =  y2[J] / (y1[J] + LAMBDAS[J]);
    temp_e_val[J]   =  exp(-LAMBDAS[J] * TAULS[J]);

    e1[J]   =  1.0 + GAMMA[J] * temp_e_val[J];  //e1                          
    e2[J]   =  1.0 - GAMMA[J] * temp_e_val[J];  //e2                      
    e3[J]   =  GAMMA[J] + temp_e_val[J];        //e3                        
    e4[J]   =  GAMMA[J] - temp_e_val[J];        //e4
  }


  J = 0;
  for(L=1; L<2*NEW_NLAYER-1; L+=2)
  {
    // HERE ARE THE EVEN MATRIX ELEMENTS
    A[L]   =  e2[J+1] * e1[J]   - e3[J] * e4[J+1];
    B[L]   =  e2[J+1] * e2[J]   - e4[J+1] * e4[J];
    D[L]   =  e1[J+1] * e4[J+1] - e3[J+1] * e2[J+1];
  
    // HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP. 
    A[L+1] =  e2[J]   * e3[J]   - e1[J]   * e4[J];
    B[L+1] =  e1[J+1] * e1[J]   - e3[J+1] * e3[J]; 
    D[L+1] =  e3[J]   * e4[J+1] - e1[J]   * e2[J+1];
    J = J + 1;

  }



  // HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
  // BEGINNING OF THE TRIDIAGONAL SOLUTION DDINITIONS. I ASSUME
  // NO DIFFUSE RADIATION IS INCIDENT AT UPPER BOUNDARY.
  A[0] = 0.0;
  B[0] = e1[0];
  D[0] = -e2[0];

  A[2*NEW_NLAYER-1] = e1[NEW_NLAYER-1] - RSFX * e3[NEW_NLAYER-1];
  B[2*NEW_NLAYER-1] = e2[NEW_NLAYER-1] - RSFX * e4[NEW_NLAYER-1];
  D[2*NEW_NLAYER-1] = 0.0;

  // This is the part of the code that solves for the blackbody stuff
  for(J=0; J<NEW_NLAYER; J++)
  {
    if(0 >= J-1)
    {
      KINDEX = 0;
    }
	  else
    {
      KINDEX = J-1;
    }
    
    temp_val_1 = (2.0 * h_constant * (NU * NU * NU)) / (CLIGHT * CLIGHT);
    temp_val_2 = exp(h_constant * NU / (bolz_constant * TEMPS[J])) - 1.0;

    B0[J] = temp_val_1 * (1.0 / temp_val_2);
    B1[J] = (B0[J] - B0[KINDEX]) / TAULS[J];
  }

  // The very top of the atmosphere is isothermal
  B1[0] = 0;

  // This solves for the C values in the toon code
  for(J=0; J<NEW_NLAYER; J++)
  {
    if(0 > J-1)
    {
      KINDEX = 0;
    }
    else
    {
      KINDEX = J-1;
    }

    temp_gamma_val[J]   = 1.0 / (y1[J] + y2[J]);

    CP[J]  = (B0[KINDEX] + B1[J] * temp_gamma_val[J]) * 2.0 * PI * mu_1;
    CPB[J] = CP[J] + B1[J] * TAULS[J] * 2.0 * PI * mu_1;

    CM[J]  = (B0[KINDEX] - B1[J] * temp_gamma_val[J]) * 2.0 * PI * mu_1;
    CMB[J] = CM[J] + B1[J] * TAULS[J] * 2.0 * PI * mu_1;
  }

  J = 0;
  for(L=1; L<2*NEW_NLAYER-1; L+=2)
  {
    E[L]   = (CP[J+1] - CPB[J]) * e2[J+1] + (CMB[J] - CM[J+1]) * e4[J+1];
    E[L+1] = e3[J] * (CP[J+1] - CPB[J]) + e1[J] * (CMB[J] - CM[J+1]);
    J = J + 1;
  }


  // HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
  // BEGINNING OF THE TRIDIAGONAL SOLUTION DDINITIONS. I ASSUME NO
  // DIFFUSE RADIATION IS INCIDENT AT THE TOP.

  E[0] = -CM[0];
  E[2*NEW_NLAYER-1]  = SFCS_HEMISPHERIC + RSFX * CMB[NEW_NLAYER-1] - CPB[NEW_NLAYER-1];
  DS[2*NEW_NLAYER-1] = E[2*NEW_NLAYER-1] / B[2*NEW_NLAYER-1];
  AS[2*NEW_NLAYER-1] = A[2*NEW_NLAYER-1] / B[2*NEW_NLAYER-1];

  //********************************************
  //*     WE SOLVE THE TRIDIAGONAL EQUATIONS   *
  //********************************************

  for(L=2; L<2*NEW_NLAYER+1; L++)
  {
    X[2*NEW_NLAYER-L]  = 1.0 / (B[2*NEW_NLAYER-L] - D[2*NEW_NLAYER-L] * AS[2*NEW_NLAYER-L+1]);
    AS[2*NEW_NLAYER-L] = A[2*NEW_NLAYER-L] * X[2*NEW_NLAYER-L];
    DS[2*NEW_NLAYER-L] = (E[2*NEW_NLAYER-L] - (D[2*NEW_NLAYER-L] * DS[2*NEW_NLAYER-L+1])) * X[2*NEW_NLAYER-L];
  }
  
  Y[0] = DS[0];
  for(L=1; L<2*NEW_NLAYER; L++)
  {
    Y[L] = DS[L] - AS[L] * Y[L-1];
  }

  //********************************************
  //*    Source Function Technique Solution    *
  //********************************************
  for(J=1; J<NEW_NLAYER+1; J++)
  {
    SOURCE_Y1[J-1] = Y[2*J-2];
    SOURCE_Y2[J-1] = Y[2*J-1];
  }


  for(J=0; J<NEW_NLAYER; J++)
  {
    if(0 > J-1)
    {
      KINDEX = 0;
    }
    else
    {
      KINDEX = J-1;
    }


    SOURCE_G[J] = (SOURCE_Y1[J] + SOURCE_Y2[J]) * (1.0/mu_1 - LAMBDAS[J]);
    SOURCE_H[J] = (SOURCE_Y1[J] - SOURCE_Y2[J]) * GAMMA[J] * (1.0/mu_1 + LAMBDAS[J]);
    SOURCE_J[J] = (SOURCE_Y1[J] + SOURCE_Y2[J]) * GAMMA[J] * (1.0/mu_1 + LAMBDAS[J]);
    SOURCE_K[J] = (SOURCE_Y1[J] - SOURCE_Y2[J]) * (1.0/mu_1 - LAMBDAS[J]);

    source_temp[J] = (1.0 / (y1[J] + y2[J])) - mu_1;

    ALPHA_1[J]     = 2.0 * PI * (B0[J] + (B1[J] * source_temp[J]));
    ALPHA_2[J]     = 2.0 * PI * (B1[J]);

    SIGMA_1[J] = 2.0 * PI * (B0[J] - (B1[J] * source_temp[J]));
    SIGMA_2[J] = 2.0 * PI * B1[J];
  }

  
  INTENSITY_DOWN[0] = BB_TOP_OF_ATM * exp(-TAULS[0]) + \
                      SOURCE_J[0]/(LAMBDAS[0] + 1.0) * (1.0 - exp(-TAULS[0] * (LAMBDAS[0]+1.0) )) + \
                      SOURCE_K[0]/(LAMBDAS[0] - 1.0) * (exp(-TAULS[0]) - exp(-TAULS[0]*LAMBDAS[0])) + \
                      SIGMA_1[0] * (1.0 - exp(-TAULS[0])) + \
                      SIGMA_2[0] * (exp(-TAULS[0]) + TAULS[0] - 1.0);


  // Do the downward intensity first
  for(J=1; J<NEW_NLAYER; J++)
  {
    INTENSITY_DOWN[J] = INTENSITY_DOWN[J-1] * exp(-TAULS[J]) + \
                    SOURCE_J[J]/(LAMBDAS[J] + 1.0) * (1.0 - exp(-TAULS[J] * (LAMBDAS[J]+1.0) )) + \
                    SOURCE_K[J]/(LAMBDAS[J] - 1.0) * (exp(-TAULS[J]) - exp(-TAULS[J]*LAMBDAS[J])) + \
                    SIGMA_1[J] * (1.0 - exp(-TAULS[J])) + \
                    SIGMA_2[J] * (exp(-TAULS[J]) + TAULS[J] - 1.0);
  }

  INTENSITY_UP[NEW_NLAYER-1] = 2.0 * BB_BOTTOM_OF_ATM * EMIS * PI;

  // Calculate the upward intensity next
  for(Z=1; Z<NEW_NLAYER; Z++)
  {
    J = NEW_NLAYER - Z - 1;

    INTENSITY_UP[J] = INTENSITY_UP[J+1] * exp(-TAULS[J+1]) + \
                      SOURCE_G[J] / (LAMBDAS[J]-1.0) * (exp(-TAULS[J+1])-exp(-TAULS[J+1] * (LAMBDAS[J]))) + \
                      SOURCE_H[J] / (LAMBDAS[J]+1.0) * (1.0 - exp(-TAULS[J+1] * (LAMBDAS[J] + 1.0))) + \
                      ALPHA_1[J] * (1.0 - exp(-TAULS[J+1])) + \
                      ALPHA_2[J] * (1.0 - ((TAULS[J+1] + 1.0) * (exp(-TAULS[J+1]))));

    HEMISPHERIC_SOURCE_FNC[J] = INTENSITY_UP[J];
  }




  //***************************************************************************************
  //***************************************************************************************
  //***************************************************************************************
  //***************************************************************************************
  //***************************************************************************************




  //**************************************************
  //*               Short Wave Solution               *
  //**************************************************

  // Boundary conditions
  mu_1 = 0.577350;
  EMIS = 1.0;
  RSFX = 1.0;
  SFCS_QUADRATURE = RSFX * mu_0 * exp(-(TAUCS[NEW_NLAYER-1]) / mu_0) * PI * FLUX_SURFACE_QUADRATURE;

  // HERE WE FIND LAYER PROPERTIES FOLLOWING GENERAL SCHEME
  // OF MEADOR AND WEAVOR. THEN WE SET UP LAYER PROPERTIES
  // NEEDED FOR MATRIX.

  for(J=0; J<NEW_NLAYER; J++)
  {
    y1[J] = 0.86602540378 * (2.0 - (W0[J] * (1.0 + G0[J])));
    y2[J] = 0.86602540378 * (W0[J] * (1.0 - G0[J]));

    y3[J] = 0.5 * (1.0 - 1.73205 * mu_0 * G0[J]);
    y4[J] = 1.0 - (y3[J]);

    LAMBDAS[J]    =  sqrt(fabs(pow(y1[J], 2.0) - pow(y2[J], 2.0)));
    GAMMA[J]  =  y2[J] / (y1[J] + LAMBDAS[J]);
    temp_e_val[J]   =  exp(-LAMBDAS[J] * TAULS[J]);

    e1[J]   =  1.0 + GAMMA[J] * temp_e_val[J];  //e1                          
    e2[J]   =  1.0 - GAMMA[J] * temp_e_val[J];  //e2                      
    e3[J]   =  GAMMA[J] + temp_e_val[J];        //e3                        
    e4[J]   =  GAMMA[J] - temp_e_val[J];        //e4
  }


  J = 0;
  for(L=1; L<2*NEW_NLAYER-1; L+=2)
  {
    // HERE ARE THE EVEN MATRIX ELEMENTS
    A[L]   =  e2[J+1] * e1[J]   - e4[J+1] * e3[J];
    B[L]   =  e2[J+1] * e2[J]   - e4[J+1] * e4[J];
    D[L]   =  e1[J+1] * e4[J+1] - e3[J+1] * e2[J+1];
  
    // HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP.
    A[L+1] =  e2[J]   * e3[J]   - e1[J]   * e4[J];
    B[L+1] =  e1[J+1] * e1[J]   - e3[J+1] * e3[J]; 
    D[L+1] =  e3[J]   * e4[J+1] - e1[J]   * e2[J+1];
    J = J + 1;
   }

  // HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
  // BEGINNING OF THE TRIDIAGONAL SOLUTION DDINITIONS. I ASSUME
  // NO DIFFUSE RADIATION IS INCIDENT AT UPPER BOUNDARY.
  A[0] = 0.0;
  B[0] = e1[0];
  D[0] = -e2[0];

  A[2*NEW_NLAYER-1] = e1[NEW_NLAYER-1] - RSFX * e3[NEW_NLAYER-1];
  B[2*NEW_NLAYER-1] = e2[NEW_NLAYER-1] - RSFX * e4[NEW_NLAYER-1];
  D[2*NEW_NLAYER-1] = 0.0;

  // This solves for the C values in the toon code
  for(J=0; J<NEW_NLAYER; J++)
  {
    if(0 > J-1)
    {
      KINDEX = 0;
    }
    else
    {
      KINDEX = J-1;
    }

    temp_gamma_val[J]   = 1.0 / (y1[J] + y2[J]);

    CP[J] = W0[J] * PI * FLUX_SURFACE_QUADRATURE * \
             exp(-(TAUCS[J]) / mu_0) * \
             (((y1[J] - 1.0 / mu_0) * y3[J]) + (y4[J] * y2[J])) / \
             (pow(LAMBDAS[J], 2.0) - (1.0 / pow(mu_0, 2.0)));


    CPB[J] = W0[J] * PI * FLUX_SURFACE_QUADRATURE * \
             exp(-(TAUCS[J] + TAULS[J]) / mu_0) *   \
             (((y1[J] - 1.0 / mu_0) * y3[J]) + (y4[J] * y2[J])) / \
             (pow(LAMBDAS[J], 2.0) - (1.0 / pow(mu_0, 2.0)));


    CM[J] = W0[J] * PI * FLUX_SURFACE_QUADRATURE * \
         exp(-(TAUCS[J]) / mu_0) *   \
         (((y1[J] + 1.0 / mu_0) * y4[J]) + (y2[J] * y3[J])) / \
         (pow(LAMBDAS[J], 2.0) - (1.0 / pow(mu_0, 2.0)));


    CMB[J] = W0[J] * PI * FLUX_SURFACE_QUADRATURE * \
             exp(-(TAUCS[J] + TAULS[J]) / mu_0) *   \
             (((y1[J] + 1.0 / mu_0) * y4[J]) + (y2[J] * y3[J])) / \
             (pow(LAMBDAS[J], 2.0) - (1.0 / pow(mu_0, 2.0)));
  }


  J = 0;
  for(L=1; L<2*NEW_NLAYER-1; L+=2)
  {
    E[L]   = (CP[J+1] - CPB[J]) * e2[J+1] + (CMB[J] - CM[J+1]) * e4[J+1];
    E[L+1] = e3[J] * (CP[J+1] - CPB[J]) + e1[J] * (CMB[J] - CM[J+1]);
    J = J + 1;
  }

  // HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
  // BEGINNING OF THE TRIDIAGONAL SOLUTION DEFINITIONS. I ASSUME NO
  // DIFFUSE RADIATION IS INCIDENT AT THE TOP.

  E[0] = -CM[0];
  E[2*NEW_NLAYER-1]  = SFCS_HEMISPHERIC + RSFX * CMB[NEW_NLAYER-1] - CPB[NEW_NLAYER-1]  ;
  DS[2*NEW_NLAYER-1] = E[2*NEW_NLAYER-1] / B[2*NEW_NLAYER-1];
  AS[2*NEW_NLAYER-1] = A[2*NEW_NLAYER-1] / B[2*NEW_NLAYER-1];

  //********************************************
  //*     WE SOLVE THE TRIDIAGONAL EQUATIONS   *
  //********************************************

  for(L=2; L<2*NEW_NLAYER+1; L++)
  {
    X[2*NEW_NLAYER-L]  = 1.0 / (B[2*NEW_NLAYER-L] - D[2*NEW_NLAYER-L] * AS[2*NEW_NLAYER-L+1]);
    AS[2*NEW_NLAYER-L] = A[2*NEW_NLAYER-L] * X[2*NEW_NLAYER-L];
    DS[2*NEW_NLAYER-L] = (E[2*NEW_NLAYER-L] - D[2*NEW_NLAYER-L] * DS[2*NEW_NLAYER-L+1]) * X[2*NEW_NLAYER-L];
  }

  Y[0] = DS[0];
  for(L=1; L<2*NEW_NLAYER; L++)
  {
    Y[L] = DS[L] - AS[L] * Y[L-1];
  }


  //***************************************************************
  //  CALCULATE LAYER CODFICIENTS, NET FLUX AND MEAN INTENSITY
  //***************************************************************

  for(J=0; J<NEW_NLAYER; J++)
  {

    TMI[J] = (1.0 / mu_1) * (Y[2*J]*(e1[J] + e3[J]) + \
             (Y[2*J+1] * (e2[J] + e4[J])) + CPB[J] + CMB[J]);

    QUADRATURE_TWO_STREAM[J] = TMI[J];
  }

  intensity_vals[0] = fabs(HEMISPHERIC_SOURCE_FNC[0]) / (2.0 * PI);
  intensity_vals[1] = fabs(QUADRATURE_TWO_STREAM[0])  / (4.0 * PI);



  //for(J=0; J<NEW_NLAYER; J++)
  //{
  //    printf("\'NEW\', %d, %le, %le, %le, \n",  J, TEMPS[J], TAUCS[J], QUADRATURE_TWO_STREAM[0]);
  //}


  //TOTAL_FLUX = fabs(HEMISPHERIC_SOURCE_FNC[0] / (2 * PI)) + (fabs(QUADRATURE_TWO_STREAM[0]) / (4.0 * PI));
  //return 0.0;
}
