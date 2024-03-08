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
 * Need to multiply by the number density to get kappa in units of m^-1
 * --------------------------------------------------------------- */

/* ------- begin ------------ ReadOpacTable.c -------------------- */
void ReadOpacTable(struct Opac *opac, const char *filename, const char *speciesName) {
    int i, j, k, fmt;
    double junk;

    FILE *f1;
    fmt = FORMAT;
    atmos.lambda = malloc(NLAMBDA * sizeof(double));

    switch (fmt) {
        case (1): {
            printf("ERROR IN CHOOSING CASE 1");
            exit(0);
            break;
        }

        case (2): {
            opac->NP = NPRESSURE;
            opac->NT = NTEMP;
            atmos.Nlambda = NLAMBDA;

            f1 = fopen(filename, "r");
            if (f1 == NULL) {
                printf("\nreadopactable.c:\nError opening file: %s -- No such file or directory\n\n", filename);
                exit(1);
            }

            for (i = 0; i < NLAMBDA; i++) {
                fscanf(f1, "%le", &atmos.lambda[i]);
            }

            fclose(f1);
            //printf("Read opacity table (format 2) from %s\n", filename);
            break;
        }

        default: {
            printf("Invalid format for opacity table\n\n");
            exit(1);
        }
    }
    printf("Read opacity table from %s\n", filename);
}

/* ------- end -------------- ReadOpacTable.c -------------------- */

/* ------- begin ------------ FreeOpacTable.c -------------------- */

void FreeOpacTable(struct Opac* opac) {
    if (opac == NULL) return; // Safety check

    // Free the temperature array
    if (opac->T != NULL) {
        free(opac->T);
        opac->T = NULL;
    }

    // Free the pressure array
    if (opac->P != NULL) {
        free(opac->P);
        opac->P = NULL;
    }

    // Free the log pressure array
    if (opac->Plog10 != NULL) {
        free(opac->Plog10);
        opac->Plog10 = NULL;
    }

    // Free the 3D kappa array
    if (opac->kappa != NULL) {
        for (int i = 0; i < NLAMBDA; i++) {
            if (opac->kappa[i] != NULL) {
                for (int j = 0; j < opac->NP; j++) {
                    free(opac->kappa[i][j]);
                }
                free(opac->kappa[i]);
            }
        }
        free(opac->kappa);
        opac->kappa = NULL;
    }

    // Free the abundance array
    if (opac->abundance != NULL) {
        for (int j = 0; j < opac->NP; j++) {
            free(opac->abundance[j]);
        }
        free(opac->abundance);
        opac->abundance = NULL;
    }
}


/* ------- end -------------- FreeOpacTable.c -------------------- */
