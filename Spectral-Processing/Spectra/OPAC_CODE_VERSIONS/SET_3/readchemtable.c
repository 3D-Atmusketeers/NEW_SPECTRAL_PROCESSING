/*----------------------- readchemtable.c ------------------------
Author: Eliza Miller-Ricci (emillerricci@cfa.harvard.edu)
Last modified: October 20, 2009
------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "input.h"
#include "opac.h"

/* --- Global variables ------------------------------------------ */

/* --- Function prototypes --------------------------------------- */

/* ---------------------------------------------------------------
 * Read in chemistry files: abundance(pressure, temperature)
 * --------------------------------------------------------------- */

/* ------- begin ------------ ReadChemTable.c -------------------- */

int ReadChemLine(FILE *f1, int i, int j, int k, double **chem_array) {
    /* Reads a chemistry line. keeps a running tally k of how many
    columns have been read in already. Also takes care of reading into
    the map_abund array. */
    fscanf(f1, "%le", &chem_array[i][j]);
    if (ABUND_TRACK_IND <= 1) {
        printf("\nMapping temperature/pressure doesn't make sense here!\n");
        exit(1);
    }
    k = k + 1;
    return k;
}

void ReadChemTable() {
    int i, j, k;
    char dum[1000];
    char speciesNames[1000];  // Buffer to store species names

    FILE *f1;

    // Open CHEM_FILE for reading
    f1 = fopen(CHEM_FILE, "r");
    if (f1 == NULL) {
        printf("\nreadchemtable.c:\nError opening file: No such file or directory\n\n");
        exit(1);
    }

    // Read species names from the first line of the file
    fgets(speciesNames, sizeof(speciesNames), f1);

    // Tokenize the species names and dynamically allocate memory for Chem structure
    char *token = strtok(speciesNames, " ");
    int speciesCount = 0;
    while (token != NULL) {
        chem.species[speciesCount++] = strdup(token);
        token = strtok(NULL, " ");
    }
    chem.numSpecies = speciesCount;

    // Allocate memory for Chem structure based on the number of species
    chem.T = malloc(NTEMP * sizeof(double));
    chem.P = malloc(NPRESSURE * sizeof(double));
    
    // Allocate memory for each species array
    chem.speciesData = malloc(chem.numSpecies * sizeof(double *));
    for (i = 0; i < chem.numSpecies; i++) {
        chem.speciesData[i] = malloc(NPRESSURE * sizeof(double *));
        for (j = 0; j < NPRESSURE; j++) {
            chem.speciesData[i][j] = malloc(NTEMP * sizeof(double));
        }
    }

    // Print the species names
    printf("Species Names:\n");
    for (i = 0; i < chem.numSpecies; i++) {
        printf("%s ", chem.species[i]);
    }
    printf("\n");

    // Read in chemistry table data
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

    // Print the first value of each species for validation
    printf("\nSpecies Data:\n");
    for (i = 0; i < chem.numSpecies; i++) {
        printf("%s: %.3e\n", chem.species[i], chem.speciesData[i][NPRESSURE - 1][NTEMP - 1]);
    }
}



/* ------- end -------------- ReadChemTable.c -------------------- */

/* ------- start ------------ FreeChemTable.c -------------------- */

/* Free allocated memory for Chem structure */
void FreeChemTable() {
    int i;
    free(chem.P);
    free(chem.T);
    // Free memory for other species arrays in Chem structure
    for (i = 0; i < NPRESSURE; i++) {
        free(chem.total[i]);
        // Free memory for other species arrays
    }
    // Free memory for species names
    for (i = 0; i < chem.numSpecies; i++) {
        free(chem.species[i]);
    }
}

/* ------- end -------------- FreeChemTable.c -------------------- */
