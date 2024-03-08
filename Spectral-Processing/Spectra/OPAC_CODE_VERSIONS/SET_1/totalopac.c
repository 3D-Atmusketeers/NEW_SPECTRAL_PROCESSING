/*----------------------- totalopac.c ----------------------------

Author: Eliza Miller-Ricci (emillerricci@cfa.harvard.edu)
Last modified: June 13, 2007

------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "input.h"
#include "include.h"
#include "constant.h"
#include "atmos.h"
#include "opac.h"
#include "nrutil.h"

/* --- Global variables ------------------------------------------ */
struct Opac** opacArray = NULL;
int opacCount = 0;

extern struct Atmos atmos;
extern struct Chem chem; // Add declaration for chem structure

/* --- Function prototypes --------------------------------------- */
void Locate(int n, double *array, double value, int *ilow);
void ReadOpacTable(struct Opac *opac, const char *filename, const char *speciesName);
double lint2D(double x1, double x2, double y1, double y2, double z1, double z2, double z3, double z4, double x, double y);
void FreeOpacTable(struct Opac* opac);
void ReadChemTable();
void FreeChemTable();

void InitializeOpac(struct Opac* opac, const char* speciesName) {
    // Error checking omitted for brevity

    // Allocate memory for T, P, and Plog10
    opac->T = malloc(NTEMP * sizeof(double));
    opac->P = malloc(NPRESSURE * sizeof(double));
    opac->Plog10 = malloc(NPRESSURE * sizeof(double));

    // Correctly allocating kappa to be a 3D array
    opac->kappa = malloc(NLAMBDA * sizeof(double**));
    for (int i = 0; i < NLAMBDA; i++) {
        opac->kappa[i] = malloc(NPRESSURE * sizeof(double*));
        for (int j = 0; j < NPRESSURE; j++) {
            opac->kappa[i][j] = malloc(NTEMP * sizeof(double));
        }
    }

    // Allocate memory for abundance and initialize it based on species data
    opac->abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
    for (int j = 0; j < NPRESSURE; j++) {
        for (int k = 0; k < NTEMP; k++) {
            // Dynamically select the right abundance data based on the species name
            for (int speciesIdx = 0; speciesIdx < chem.numSpecies; speciesIdx++) {
                if (strcmp(speciesName, chem.species[speciesIdx]) == 0) {
                    opac->abundance[j][k] = chem.speciesData[speciesIdx][j][k];
                    break; // Break the loop once the species data is found
                }
            }
        }
    }
    
    // Set the species name
    strncpy(opac->speciesName, speciesName, MAX_NAME_LENGTH - 1);
    opac->speciesName[MAX_NAME_LENGTH - 1] = '\0'; // Ensure null termination

    // Read the opacity table which will populate kappa (and possibly adjust based on abundance)
    char filename[256];
    snprintf(filename, sizeof(filename), "DATA/SET_1/opac%s.dat", speciesName);
    ReadOpacTable(opac, filename, speciesName); // Call the function without 'struct' keyword

    printf("Initialized and added species '%s'.\n", speciesName);

    // Assuming opacArray and opacCount are globally accessible
    opacArray[opacCount++] = opac;
}

/* ---------------------------------------------------------------
 * Computes the total opacity due to all of the atmospheric 
 * constituents.
 * --------------------------------------------------------------- */

/* ------- begin ------------ TotalOpac.c ------------------------ */
void TotalOpac() {
  double **opac_CIA_H2H2, **opac_CIA_H2He, **opac_CIA_H2H, **opac_CIA_H2CH4, **opac_CIA_CH4Ar,
         **opac_CIA_CH4CH4, **opac_CIA_CO2CO2, **opac_CIA_HeH, **opac_CIA_N2CH4, **opac_CIA_N2H2,
         **opac_CIA_N2N2, **opac_CIA_O2CO2, **opac_CIA_O2N2, **opac_CIA_O2O2, **opac_CIA_Hel;
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

  /* Read Chemistry Table */
  ReadChemTable();
  printf("ReadChemTable done\n");

  exit(0);


  // Structure to hold species names
  typedef struct {
      char names[MAX_OPACITY_SPECIES][MAX_NAME_LENGTH];
      int count;
  } SpeciesList;

  FILE* file = fopen("input.h", "r");
  if (!file) {
      perror("Failed to open file");
      exit(0);
  }
  printf("File opened successfully.\n");

  SpeciesList speciesList = {.count = 0};
  char line[256];
  char* startPattern = "#define ";
  char* endPattern = ".dat\"";

  while (fgets(line, sizeof(line), file)) {
      //printf("Reading line: %s", line); // Debug: print the current line
      char* startPos = strstr(line, startPattern);
      char* endPos = strstr(line, endPattern);
      if (startPos && endPos) {
          //printf("Found a line with species data.\n"); // Debug: found a line of interest
          startPos += strlen(startPattern); // Move past "#define "
          // Calculate name length, ensuring we don't exceed buffer sizes
          ptrdiff_t nameLength = endPos - startPos + strlen(".dat"); // Include ".dat" in the name
          if (nameLength > 0 && nameLength < MAX_NAME_LENGTH + 4) { // "+4" to account for "opac" and ".dat"
              char tempName[MAX_NAME_LENGTH + 20]; // Temporary buffer to hold the full filename
              strncpy(tempName, startPos, nameLength);
              tempName[nameLength] = '\0'; // Ensure null-termination
              
              if (strstr(tempName, "chem") == NULL && strstr(tempName, "CIA") == NULL) {
                  char* lastSlash = strrchr(tempName, '/'); // Find the last slash to isolate the filename
                  if (lastSlash != NULL) {
                      char* speciesName = lastSlash + 5; // Skip past the slash and "opac"
                      int speciesNameLength = strlen(speciesName) - 4; // Remove the ".dat" part
                      if (speciesNameLength > 0 && speciesNameLength < MAX_NAME_LENGTH) {
                          strncpy(speciesList.names[speciesList.count], speciesName, speciesNameLength);
                          speciesList.names[speciesList.count][speciesNameLength] = '\0'; // Ensure null-termination
                          //printf("Extracted species name: %s\n", speciesList.names[speciesList.count]); // Debug
                          speciesList.count++;
                          if (speciesList.count >= MAX_OPACITY_SPECIES) {
                              printf("Reached maximum species count.\n"); // Debug
                              break; // Avoid exceeding array bounds
                          }
                      }
                  }
              } else {
                  //printf("Skipped species name due to filter: %s\n", tempName); // Debug
              }
          } else {
              printf("Name length issue: Length is %td, expected < %d.\n", nameLength, MAX_NAME_LENGTH); // Debug
          }
      }
  }

  fclose(file);
  //printf("File closed.\n");

  // Print the species names
  printf("Extracted Species Names:\n");
  for (int i = 0; i < speciesList.count; ++i) {
      printf("%d: %s\n", i + 1, speciesList.names[i]);
  }

  for (int i = 0; i < speciesList.count; ++i) {
      // Check if opacCount has not exceeded the limit
      if (opacCount >= MAX_OPACITY_SPECIES) {
          fprintf(stderr, "Maximum number of species exceeded.\n");
          return;
      }

      // Allocate memory for a new Opac structure
      struct Opac* newOpac = (struct Opac*)malloc(sizeof(struct Opac));
      if (!newOpac) {
          fprintf(stderr, "Memory allocation for new Opac failed.\n");
          return;
      }

      // Initialize the new Opac structure for the current species
      InitializeOpac(newOpac, speciesList.names[i]); // Make sure InitializeOpac function is implemented correctly

      // Add the newly initialized Opac to the global array
      opacArray[opacCount++] = newOpac;
  }



  /*
  // Fill in collision-induced opacities 
  f1 = fopen("DATA/SET_1/opacCIA_low_temp.dat", "r");
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


      opacCIA.kappa[i][j][k] += opac_CIA_Hel[k][i] *
	       chem.H[j][k] * chem.P[j] / (KBOLTZMANN * chem.T[k]) *  
	       chem.el[j][k] * chem.P[j] / (KBOLTZMANN * chem.T[k]);
      }
    }
  }
  */

  fclose(f1);
	
  /* Rayleigh scattering */
  /* (Polarizabilities from the CRC Handbook) */

  for (i = 0; i < NLAMBDA; i++) {
      for (j = 0; j < NPRESSURE; j++) {
          for (k = 0; k < NTEMP; k++) {
              double rayleigh_contrib = 0.0;

              // Check and accumulate contributions from each species
              if (chem.speciesData != NULL) {
                  int speciesIdx;
                  for (speciesIdx = 0; speciesIdx < chem.numSpecies; speciesIdx++) {
                      if (strcmp(chem.species[speciesIdx], "H2") == 0) {
                          if (chem.speciesData[speciesIdx] != NULL) {
                              rayleigh_contrib += (8.0 * PI / 3.0) * SQ(0.80e-30) *
                                                  SQ(2.0 * PI / atmos.lambda[i]) * SQ(2.0 * PI / atmos.lambda[i]) *
                                                  chem.speciesData[speciesIdx][j][k] * chem.P[j] / (KBOLTZMANN * chem.T[k]);
                          } else {
                              printf("Error: Species H2 data is missing.\n");
                              exit(1);
                          }
                      } else if (strcmp(chem.species[speciesIdx], "He") == 0) {
                          if (chem.speciesData[speciesIdx] != NULL) {
                              rayleigh_contrib += (8.0 * PI / 3.0) * SQ(0.21e-30) *
                                                  SQ(2.0 * PI / atmos.lambda[i]) * SQ(2.0 * PI / atmos.lambda[i]) *
                                                  chem.speciesData[speciesIdx][j][k] * chem.P[j] / (KBOLTZMANN * chem.T[k]);
                          } else {
                              printf("Error: Species He data is missing.\n");
                              exit(1);
                          }
                      } else if (strcmp(chem.species[speciesIdx], "H2O") == 0) {
                          if (chem.speciesData[speciesIdx] != NULL) {
                              rayleigh_contrib += (8.0 * PI / 3.0) * SQ(1.45e-30) *
                                                  SQ(2.0 * PI / atmos.lambda[i]) * SQ(2.0 * PI / atmos.lambda[i]) *
                                                  chem.speciesData[speciesIdx][j][k] * chem.P[j] / (KBOLTZMANN * chem.T[k]);
                          } else {
                              printf("Error: Species H2O data is missing.\n");
                              exit(1);
                          }
                      } else if (strcmp(chem.species[speciesIdx], "CO") == 0) {
                          if (chem.speciesData[speciesIdx] != NULL) {
                              rayleigh_contrib += (8.0 * PI / 3.0) * SQ(1.95e-30) *
                                                  SQ(2.0 * PI / atmos.lambda[i]) * SQ(2.0 * PI / atmos.lambda[i]) *
                                                  chem.speciesData[speciesIdx][j][k] * chem.P[j] / (KBOLTZMANN * chem.T[k]);
                          } else {
                              printf("Error: Species CO data is missing.\n");
                              exit(1);
                          }
                      } else if (strcmp(chem.species[speciesIdx], "CO2") == 0) {
                          if (chem.speciesData[speciesIdx] != NULL) {
                              rayleigh_contrib += (8.0 * PI / 3.0) * SQ(2.91e-30) *
                                                  SQ(2.0 * PI / atmos.lambda[i]) * SQ(2.0 * PI / atmos.lambda[i]) *
                                                  chem.speciesData[speciesIdx][j][k] * chem.P[j] / (KBOLTZMANN * chem.T[k]);
                          } else {
                              printf("Error: Species CO2 data is missing.\n");
                              exit(1);
                          }
                      } else if (strcmp(chem.species[speciesIdx], "NH3") == 0) {
                          if (chem.speciesData[speciesIdx] != NULL) {
                              rayleigh_contrib += (8.0 * PI / 3.0) * SQ(2.26e-30) *
                                                  SQ(2.0 * PI / atmos.lambda[i]) * SQ(2.0 * PI / atmos.lambda[i]) *
                                                  chem.speciesData[speciesIdx][j][k] * chem.P[j] / (KBOLTZMANN * chem.T[k]);
                          } else {
                              printf("Error: Species NH3 data is missing.\n");
                              exit(1);
                          }
                      } else if (strcmp(chem.species[speciesIdx], "CH4") == 0) {
                          if (chem.speciesData[speciesIdx] != NULL) {
                              rayleigh_contrib += (8.0 * PI / 3.0) * SQ(2.59e-30) *
                                                  SQ(2.0 * PI / atmos.lambda[i]) * SQ(2.0 * PI / atmos.lambda[i]) *
                                                  chem.speciesData[speciesIdx][j][k] * chem.P[j] / (KBOLTZMANN * chem.T[k]);
                          } else {
                              printf("Error: Species CH4 data is missing.\n");
                              exit(1);
                          }
                      }
                      // Add similar checks for other species as necessary...
                  }
              } else {
                  printf("Error: Species data is missing.\n");
                  exit(1);
              }

              // Accumulate contribution to opacity
              opac->kappa[i][j][k] += rayleigh_contrib;
          }
      }
  }
  
  /* Free uneeded opacity structures and chemistry table */

  //FreeOpacTable(opacC2H2);
  //FreeOpacTable(opacCH4);
  //FreeOpacTable(opacCO);
  //FreeOpacTable(opacCO2);
  //FreeOpacTable(opacFeH);
  //FreeOpacTable(opacH2O);
  //FreeOpacTable(opacH2S);
  //FreeOpacTable(opacHCN);
  //FreeOpacTable(opacK);
  //FreeOpacTable(opacNa);
  //FreeOpacTable(opacNH3);
  //FreeOpacTable(opacTiO);
  //FreeOpacTable(opacVO);

  //FreeOpacTable(opacCIA);
  //FreeOpacTable(opacscat);
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



