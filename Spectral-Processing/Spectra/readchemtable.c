#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "input.h"
#include "opac.h"

/* --- Global variables ------------------------------------------ */
struct Chem chem;

/* --- Function prototypes --------------------------------------- */
int ReadChemLine(FILE *f1, int i, int j, int k, double **chem_array);

/* ---------------------------------------------------------------
 * Read in chemistry files: abundance(pressure, temperature)
 * --------------------------------------------------------------- */

/* ------- begin ------------ ReadChemTable.c -------------------- */

void ReadChemTable() {
    int i, j, k;
    char speciesNames[1000];

    FILE *f1 = fopen(CHEM_FILE, "r");
    if (f1 == NULL) {
        fprintf(stderr, "Error opening file: %s\n", CHEM_FILE);
        exit(EXIT_FAILURE);
    }

    fgets(speciesNames, sizeof(speciesNames), f1);

    char *token = strtok(speciesNames, " ");
    int speciesCount = 0;
    while (token != NULL) {
        chem.species[speciesCount++] = strdup(token);
        token = strtok(NULL, " ");
    }
    chem.numSpecies = speciesCount;

    chem.T = malloc(NTEMP * sizeof(double));
    chem.P = malloc(NPRESSURE * sizeof(double));
    chem.species = malloc(chem.numSpecies * sizeof(char *));
    chem.speciesData = malloc(chem.numSpecies * sizeof(double *));
    for (i = 0; i < chem.numSpecies; i++) {
        chem.speciesData[i] = malloc(NPRESSURE * sizeof(double *));
        for (j = 0; j < NPRESSURE; j++) {
            chem.speciesData[i][j] = malloc(NTEMP * sizeof(double));
        }
    }

    for (i = NPRESSURE - 1; i >= 0; i--) {
        fscanf(f1, "%le", &chem.P[i]);
        for (j = NTEMP - 1; j >= 0; j--) {
            fscanf(f1, "%le", &chem.T[j]);
            k = 0;
            for (int speciesIdx = 0; speciesIdx < chem.numSpecies; speciesIdx++) {
                k = ReadChemLine(f1, i, j, k, chem.speciesData[speciesIdx]);
            }
        }
    }
    fclose(f1);

    printf("Species Names:\n");
    for (i = 0; i < chem.numSpecies; i++) {
        printf("%s ", chem.species[i]);
    }
    printf("\n");

    printf("\nSpecies Data:\n");
    for (i = 0; i < chem.numSpecies; i++) {
        printf("%s: %.3e\n", chem.species[i], chem.speciesData[i][NPRESSURE - 1][NTEMP - 1]);
    }
}

int ReadChemLine(FILE *f1, int i, int j, int k, double **chem_array) {
    fscanf(f1, "%le", &chem_array[i][j]);
    if (ABUND_TRACK_IND <= 1) {
        fprintf(stderr, "Mapping temperature/pressure doesn't make sense here!\n");
        exit(EXIT_FAILURE);
    }
    k = k + 1;
    return k;
}

/* ------- end -------------- ReadChemTable.c -------------------- */

/* ------- start ------------ FreeChemTable.c -------------------- */

void FreeChemTable() {
    free(chem.P);
    free(chem.T);
    for (int i = 0; i < chem.numSpecies; i++) {
        free(chem.species[i]);
        for (int j = 0; j < NPRESSURE; j++) {
            free(chem.speciesData[i][j]);
        }
        free(chem.speciesData[i]);
    }
    free(chem.species);
    free(chem.speciesData);
}

/* ------- end -------------- FreeChemTable.c -------------------- */
