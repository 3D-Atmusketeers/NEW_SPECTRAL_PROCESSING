#include <stdio.h>


double Planck(double T, double lambda);

void two_stream(int num_tau_layers, int NLAYER, int kmin, double *w0_array, double *g0_array, \
                 const double *temperature_array, const double *tau_array, \
                 double NU, double NU_BIN, double incident_frac, double *dtau_array, \
                 double intensity_vals[])
{
  double mu_1;  // Param for Quadrature or Hemispheric Mean constant
  double mu_0;  // Incident angle of the solar beam
  double mu;    // Never used. The angle of the upwards intensity

  // These are indexing values
  int J, L, KINDEX, Z, M, g;

  // These are just constants
  double bolz_constant = 1.380649e-23;
  double h_constant    = 6.62607015e-34;

  // These are boundary conditions
  // They vary for long wave and short wave
  double EMIS;
  double RSFX;
  double SFCS_HEMISPHERIC;

  // The direct radiation from the star
  double DIRECT_HEMISPHERIC[num_tau_layers - kmin + 1];
  double DIRECT_QUADRATURE[num_tau_layers - kmin + 1];

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
  double W0[num_tau_layers - kmin + 1];
  double G0[num_tau_layers - kmin + 1];

  // Scattering and atmosphere parameters
  double TEMPS[num_tau_layers - kmin + 2];

  // Set the optical depths
  double TAUCS[num_tau_layers - kmin + 2];
  double TAULS[num_tau_layers - kmin + 1];

  // How to generate the matrix from the toon paper
  double y1[num_tau_layers - kmin + 1];
  double y2[num_tau_layers - kmin + 1];
  double y3[num_tau_layers - kmin + 1];
  double y4[num_tau_layers - kmin + 1];

  // More Matrix stuff
  double LAMBDAS[num_tau_layers - kmin + 1];
  double GAMMA[num_tau_layers - kmin + 1];
  double temp_e_val[num_tau_layers - kmin + 1];
  double exptrm_positive[num_tau_layers - kmin + 1];
  double exptrm_negative[num_tau_layers - kmin + 1];
  double e1[num_tau_layers - kmin + 1];
  double e2[num_tau_layers - kmin + 1];
  double e3[num_tau_layers - kmin + 1];
  double e4[num_tau_layers - kmin + 1];

  // Matrix Coefficients
  double A[2 * (num_tau_layers - kmin + 1)];
  double B[2 * (num_tau_layers - kmin + 1)];
  double D[2 * (num_tau_layers - kmin + 1)];
  double E[2 * (num_tau_layers - kmin + 1)];

  // Solution Coefficients
  double AS[2 * (num_tau_layers - kmin + 1)];
  double DS[2 * (num_tau_layers - kmin + 1)];
  double X[2 * (num_tau_layers - kmin + 1)];
  double Y[2 * (num_tau_layers - kmin + 1)];

  // Temporary values, makes math easier for me
  double temp_gamma_val[num_tau_layers - kmin + 1];

  // Planck Intensities
  // B1 is the slope of intensities between layers
  double B0[num_tau_layers - kmin + 2];
  double B1[num_tau_layers - kmin + 1];
  double TEMP[num_tau_layers - kmin + 1];


  // Stuff for solving the blackbody equation
  double Bnu, twohnu3_c2, hc_Tkla;
  double temp_val_1, temp_val_2;

  // Matrix values, related to the intensity at each layer
  // These set differently for the quadrature and hemipsheric solutions
  double CP[num_tau_layers - kmin + 1];
  double CPB[num_tau_layers - kmin + 1];
  double CM[num_tau_layers - kmin + 1];
  double CMB[num_tau_layers - kmin + 1];

  // The Source Function Variables
  double SOURCE_G[num_tau_layers - kmin + 1];
  double SOURCE_H[num_tau_layers - kmin + 1];
  double SOURCE_J[num_tau_layers - kmin + 1];
  double SOURCE_K[num_tau_layers - kmin + 1];
  double ALPHA_1[num_tau_layers - kmin + 1];
  double ALPHA_2[num_tau_layers - kmin + 1];
  double SIGMA_1[num_tau_layers - kmin + 1];
  double SIGMA_2[num_tau_layers - kmin + 1];
  double SOURCE_Y1[num_tau_layers - kmin + 1];
  double SOURCE_Y2[num_tau_layers - kmin + 1];
  double source_temp[num_tau_layers - kmin + 1];

  double b_top;

  double gangle[] = {0.0985350858, 0.3045357266, 0.5620251898, 0.8019865821, 0.9601901429};
  double gweight[] = {0.0157479145, 0.0739088701, 0.1463869871, 0.1671746381, 0.0967815902};
  int num_gangle = 5;

  //double gangle[] = {1.0};
  //double gweight[] = {0.5};
  //int num_gangle = 1;

  // Upward and downwards intensities
  double INTENSITY_DOWN[num_tau_layers - kmin + 2][num_gangle];
  double INTENSITY_UP[num_tau_layers - kmin + 2][num_gangle];
  double TMI[num_tau_layers - kmin + 1];
  double FLUX_UP[num_tau_layers - kmin + 1];

  // Final values for long wave radiation
  double HEMISPHERIC_TWO_STREAM[num_tau_layers - kmin + 1];
  double HEMISPHERIC_SOURCE_FNC[num_tau_layers - kmin + 1];

  double TOTAL_FLUX;

  // Final values for short wave radiation
  double QUADRATURE_TWO_STREAM[num_tau_layers - kmin + 1];
  double QUADRATURE_SOURCE_FNC[num_tau_layers - kmin + 1];

  // Top layer values short wave
  double QUADRATURE_TWO_STREAM_TOP;
  double QUADRATURE_SOURCE_FNC_TOP;

  // Calculate the flux at the top of the atmosphere from the star
  temp_val_1 = (2.0 * h_constant * (NU * NU * NU)) / (CLIGHT * CLIGHT);
  temp_val_2 = exp(h_constant * NU / (bolz_constant * STELLAR_TEMP)) - 1.0;
  STELLAR_BB = temp_val_1 * (1.0 / temp_val_2);

  // Get the flux at the surface of the atmosphere
  FLUX_SURFACE_HEMISPHERIC = 0;
  FLUX_SURFACE_QUADRATURE = incident_frac * PI * STELLAR_BB * pow((R_STAR / ORB_SEP), 2.0);

  mu_0 = 1.0;

  // Assign arrays for the characteristics of the atmosphere based on the input arrays
  TAUCS[0] = 0.0;
  for (J=0; J<NLAYER; J++)
  {
    W0[J] = w0_array[J+kmin-1];
    G0[J] = g0_array[J+kmin-1];

    TAULS[J]   = dtau_array[J+kmin-1];
    TAUCS[J+1] = TAUCS[J]+dtau_array[J+kmin-1];

    TEMPS[J] = temperature_array[J+kmin-1];

    DIRECT_QUADRATURE[J]  = mu_0 * PI * FLUX_SURFACE_QUADRATURE * exp(-1.0 * (TAUCS[J] + TAULS[J]) / mu_0);
    DIRECT_HEMISPHERIC[J] = 0.0;
  }
  TEMPS[NLAYER] = TEMPS[NLAYER - 1];


  // Calculate the intensity at the top of the atmosphere
  temp_val_1 = (2.0 * h_constant * (NU * NU * NU)) / (CLIGHT * CLIGHT);
  temp_val_2 = exp(h_constant * NU / (bolz_constant * TEMPS[0])) - 1.0;
  BB_TOP_OF_ATM = temp_val_1 * (1.0 / temp_val_2);

  // Calculate the intensity at the bottom of the atmosphere
  temp_val_1 = (2.0 * h_constant * (NU * NU * NU)) / (CLIGHT * CLIGHT);
  temp_val_2 = exp(h_constant * NU / (bolz_constant * TEMPS[NLAYER])) - 1.0;
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
  RSFX = 0.0;

  SFCS_HEMISPHERIC = EMIS * PI * BB_BOTTOM_OF_ATM;

  // HERE WE FIND LAYER PROPERTIES FOLLOWING GENERAL SCHEME
  // OF MEADOR AND WEAVOR. THEN WE SET UP LAYER PROPERTIES
  // NEEDED FOR MATRIX.

  for(J=0; J<NLAYER; J++)
  {
    y1[J]    =  (2.0 - (W0[J] * (1.0 + G0[J])));
    y2[J]    =  (W0[J] * (1.0 - G0[J]));

    LAMBDAS[J]    =  sqrt(fabs(pow(y1[J], 2.0) - pow(y2[J], 2.0)));
    GAMMA[J]  =  y2[J] / (y1[J] + LAMBDAS[J]);

    temp_e_val[J]   =  exp(-LAMBDAS[J] * TAULS[J]);
    exptrm_positive[J] = exp(fmin(LAMBDAS[J] * TAULS[J], 35));
    exptrm_negative[J] = 1.0 / exptrm_positive[J];

    e1[J]   =  exptrm_positive[J] + GAMMA[J] * temp_e_val[J];  //e1
    e2[J]   =  exptrm_positive[J] - GAMMA[J] * temp_e_val[J];  //e2
    e3[J]   =  GAMMA[J]*exptrm_positive[J] + temp_e_val[J];        //e3
    e4[J]   =  GAMMA[J]*exptrm_positive[J] - temp_e_val[J];        //e4
  }


  // This is the part of the code that solves for the blackbody stuff
  for(J=0; J<NLAYER + 1; J++)
  {
    temp_val_1 = (2.0 * h_constant * (NU * NU * NU)) / (CLIGHT * CLIGHT);
    temp_val_2 = exp(h_constant * NU / (bolz_constant * TEMPS[J])) - 1.0;
    B0[J] = temp_val_1 * (1.0 / temp_val_2);
  }


  // This is the part of the code that solves for the blackbody stuff
  for(J=0; J<NLAYER; J++)
  {
    temp_val_1 = (2.0 * h_constant * (NU * NU * NU)) / (CLIGHT * CLIGHT);
    temp_val_2 = exp(h_constant * NU / (bolz_constant * TEMPS[J])) - 1.0;
    B1[J] = (B0[J+1] - B0[J]) / TAULS[J];
  }

  b_top = (1.0 - exp(-TAULS[0]/ mu_1)) * B0[0] * PI;

  // The very top of the atmosphere is isothermal
  B1[0] = (B0[1] - B0[0]) / TAULS[0];

  // This solves for the C values in the toon code
  for(J=0; J<NLAYER; J++)
  {
    temp_gamma_val[J]   = 1.0 / (y1[J] + y2[J]);

    CP[J]  = (B0[J] + B1[J] * temp_gamma_val[J]) * 2.0 * PI * mu_1;
    CM[J]  = (B0[J] - B1[J] * temp_gamma_val[J]) * 2.0 * PI * mu_1;

    CPB[J] = CP[J] + B1[J] * TAULS[J] * 2.0 * PI * mu_1;
    CMB[J] = CM[J] + B1[J] * TAULS[J] * 2.0 * PI * mu_1;
  }



  // HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
  // BEGINNING OF THE TRIDIAGONAL SOLUTION DDINITIONS. I ASSUME
  // NO DIFFUSE RADIATION IS INCIDENT AT UPPER BOUNDARY.
  A[0] = 0;
  B[0] = GAMMA[0] + 1;
  D[0] = GAMMA[0] - 1;

  A[2*NLAYER-1] = e1[NLAYER-1] - RSFX * e3[NLAYER-1];
  B[2*NLAYER-1] = e2[NLAYER-1] - RSFX * e4[NLAYER-1];
  D[2*NLAYER-1] = 0;

  for(J=0; J<NLAYER-1; J++)
  {
    L = 2*(J+1);
    // HERE ARE THE ODD MATRIX ELEMENTS
    A[L] = 2.0 * (1.0- pow(GAMMA[J], 2.0));
    B[L] = (e1[J]-e3[J]) * (GAMMA[J+1] + 1.0);
    D[L] = (e1[J]+e3[J]) * (GAMMA[J+1] - 1.0);
  }

  for(J=0; J<NLAYER-1; J++)
  {
    L = 2*J+1;
    // HERE ARE THE EVEN MATRIX ELEMENTS
    A[L] = (e1[J]+e3[J]) * (GAMMA[J+1]-1.0);
    B[L] = (e2[J]+e4[J]) * (GAMMA[J+1]-1.0);
    D[L] = 2.0 * (1.0- pow(GAMMA[J+1], 2.0));
  }


  for(J=0; J<NLAYER-1; J++)
  {
    // ODD TERMS
    L = 2*(J+1);
    E[L]   = e3[J] * (CP[J+1] - CPB[J]) + e1[J] * (CMB[J] - CM[J+1]);
  }


  for(J=0; J<NLAYER-1; J++)
  {
    // EVEN TERMS
    L = 2*J+1;
    E[L] = (GAMMA[J+1]-1.0)*(CP[J+1] - CPB[J]) + (1.0-GAMMA[J+1])*(CMB[J] - CM[J+1]);
  }

  // HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
  // BEGINNING OF THE TRIDIAGONAL SOLUTION DDINITIONS. I ASSUME NO
  // DIFFUSE RADIATION IS INCIDENT AT THE TOP.

  E[0] = b_top - CM[0];
  E[2*NLAYER-1]  = 0;

  AS[2*NLAYER-1] = A[2*NLAYER-1] / B[2*NLAYER-1];
  DS[2*NLAYER-1] = E[2*NLAYER-1] / B[2*NLAYER-1];

  //********************************************
  //*     WE SOLVE THE TRIDIAGONAL EQUATIONS   *
  //********************************************

  for(L=2; L<2*NLAYER + 1; L++)
  {
    X[2*NLAYER-L]  = 1.0 / (B[2*NLAYER-L] - D[2*NLAYER-L] * AS[2*NLAYER-L+1]);
    AS[2*NLAYER-L] = A[2*NLAYER-L] * X[2*NLAYER-L];
    DS[2*NLAYER-L] = (E[2*NLAYER-L] - (D[2*NLAYER-L] * DS[2*NLAYER-L+1])) * X[2*NLAYER-L];
  }

  Y[0] = DS[0];
  for(L=1; L<2*NLAYER; L++)
  {
    Y[L] = DS[L] - AS[L] * Y[L-1];
  }

  //********************************************
  //*    Source Function Technique Solution    *
  //********************************************
  for(J=0; J<NLAYER; J++)
  {
    SOURCE_Y1[J] = Y[2*J];
    SOURCE_Y2[J] = Y[2*J+1];
  }


  for(J=0; J<NLAYER; J++)
  {
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



  for(g=0; g<num_gangle; g++)
  {
    mu = gangle[g];
    INTENSITY_UP[NLAYER][g] = (2.0 * PI) * (B1[NLAYER-1]  + B0[NLAYER] * mu);
    INTENSITY_DOWN[0][g] =  (1. - exp(-TAULS[0] / mu)) * BB_TOP_OF_ATM * 2. * PI;



    for(J=0; J<NLAYER; J++)
    {
      INTENSITY_DOWN[J+1][g] = INTENSITY_DOWN[J][g] * exp(-TAULS[J]/mu) + \
                      SOURCE_J[J]/(LAMBDAS[J] * mu + 1.0) * (exptrm_positive[J] - exp(-TAULS[J]/mu)) + \
                      SOURCE_K[J]/(LAMBDAS[J] * mu - 1.0) * (exp(-TAULS[J]/mu) - exptrm_negative[J]) + \
                      SIGMA_1[J] * (1.0 - exp(-TAULS[J]/mu)) + \
                      SIGMA_2[J] * (mu * exp(-TAULS[J]/mu) + TAULS[J] - mu);



      Z = NLAYER - 1 - J;
      INTENSITY_UP[Z][g] = INTENSITY_UP[Z+1][g] * exp(-TAULS[Z]/mu) + \
                        SOURCE_G[Z]/(LAMBDAS[Z]*mu-1.0) * (exptrm_positive[Z]*exp(-TAULS[Z]/mu) - 1.0) + \
                        SOURCE_H[Z]/(LAMBDAS[Z]*mu+1.0) * (1.0 - exptrm_negative[Z]*exp(-TAULS[Z]/mu)) + \
                        ALPHA_1[Z] * (1.0 - exp(-TAULS[Z]/mu)) + \
                        ALPHA_2[Z] * (mu - ((TAULS[Z] + mu) * (exp(-TAULS[Z]/mu))));
      
    }
  }


  // Calculate the weighted sum for each layer
  for(int J = 0; J < NLAYER; J++)
  {
    double weighted_sum = 0.0;
    for(int g = 0; g < num_gangle; g++)
    {
      weighted_sum += INTENSITY_UP[J][g] * gweight[g];
    }
    HEMISPHERIC_SOURCE_FNC[J] = weighted_sum;
  }

  for(int J=0; J<NLAYER; J++)
  {
    printf("%d, %0.5e, \n", J, HEMISPHERIC_SOURCE_FNC[J] * 2);
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
  RSFX = 0.0;

  // HERE WE FIND LAYER PROPERTIES FOLLOWING GENERAL SCHEME
  // OF MEADOR AND WEAVOR. THEN WE SET UP LAYER PROPERTIES
  // NEEDED FOR MATRIX.

  for(J=0; J<NLAYER; J++)
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




  // HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
  // BEGINNING OF THE TRIDIAGONAL SOLUTION DDINITIONS. I ASSUME
  // NO DIFFUSE RADIATION IS INCIDENT AT UPPER BOUNDARY.
  A[0] = 0;
  B[0] = e1[0];
  D[0] = -1 * e2[0];

  A[2*NLAYER-1] = e1[NLAYER-1] - RSFX * e3[NLAYER-1];
  B[2*NLAYER-1] = e2[NLAYER-1] - RSFX * e4[NLAYER-1];
  D[2*NLAYER-1] = 0;

  for(J=0; J<NLAYER-1; J++)
  {
    L = 2*(J+1);
    // HERE ARE THE ODD MATRIX ELEMENTS
    A[L] =  e2[J] * e3[J]   - e1[J]   * e4[J];
    B[L] =  e1[J] * e1[J+1] - e3[J+1] * e3[J];
    D[L] =  e3[J] * e4[J+1] - e1[J]   * e2[J+1];
  }

  for(J=0; J<NLAYER-1; J++)
  {
    L = 2*J+1;
    // HERE ARE THE EVEN MATRIX ELEMENTS
    A[L]   =  e2[J+1] * e1[J]   - e3[J] * e4[J+1];
    B[L]   =  e2[J] * e2[J+1]   - e4[J+1] * e4[J];
    D[L]   =  e1[J+1] * e4[J+1] - e3[J+1] * e2[J+1];
  }



  // This solves for the C values in the toon code
  for(J=0; J<NLAYER; J++)
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

  for(J=0; J<NLAYER-1; J++)
  {
    L = 2*(J+1);
    E[L]   = e3[J] * (CP[J+1] - CPB[J]) + e1[J] * (CMB[J] - CM[J+1]);
  }


  for(J=0; J<NLAYER-1; J++)
  {
    L = 2*J+1;
    E[L]   = e2[J+1] * (CP[J+1] - CPB[J]) + (CMB[J] - CM[J+1]) * e4[J+1];
  }


  E[0] = -CM[0];
  E[2*NLAYER-1]  = 0;

  // HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
  // BEGINNING OF THE TRIDIAGONAL SOLUTION DEFINITIONS. I ASSUME NO
  // DIFFUSE RADIATION IS INCIDENT AT THE TOP.

  DS[2*NLAYER-1] = E[2*NLAYER-1] / B[2*NLAYER-1];
  AS[2*NLAYER-1] = A[2*NLAYER-1] / B[2*NLAYER-1];

  //********************************************
  //*     WE SOLVE THE TRIDIAGONAL EQUATIONS   *
  //********************************************

  for(L=2; L<2*NLAYER+1; L++)
  {
    X[2*NLAYER-L]  = 1.0 / (B[2*NLAYER-L] - D[2*NLAYER-L] * AS[2*NLAYER-L+1]);
    AS[2*NLAYER-L] = A[2*NLAYER-L] * X[2*NLAYER-L];
    DS[2*NLAYER-L] = (E[2*NLAYER-L] - D[2*NLAYER-L] * DS[2*NLAYER-L+1]) * X[2*NLAYER-L];
  }

  Y[0] = DS[0];
  for(L=1; L<2*NLAYER; L++)
  {
    Y[L] = DS[L] - AS[L] * Y[L-1];
  }


  //***************************************************************
  //  CALCULATE LAYER CODFICIENTS, NET FLUX AND MEAN INTENSITY
  //***************************************************************

  for(J=0; J<NLAYER; J++)
  {
    TMI[J] = ((1.0 / mu_1) * (Y[2*J]*(e1[J] + e3[J]) + \
             (Y[2*J+1] * (e2[J] + e4[J])) + CPB[J] + CMB[J]) + DIRECT_QUADRATURE[J]/mu_0);

    FLUX_UP[J] = Y[2*J]   * (exp(-LAMBDAS[J]*(TAULS[J]-0)) + (GAMMA[J]*exptrm_negative[J])) + \
                 Y[2*J+1] * (exp(-LAMBDAS[J]*(TAULS[J]-0)) - (GAMMA[J]*exptrm_negative[J])) + \
                 CPB[J];
                 
    QUADRATURE_TWO_STREAM[J] = FLUX_UP[J];
  }

  intensity_vals[0] = fabs(HEMISPHERIC_SOURCE_FNC[0]) / (PI);
  intensity_vals[1] = fabs(QUADRATURE_TWO_STREAM[0])  / (4.0 * PI);
}
