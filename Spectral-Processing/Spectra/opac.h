/* ------- file: -------------------------- opac.h ------------------ */

#ifndef __OPAC_H__
#define __OPAC_H__

#include <stdlib.h> // For size_t definition

// Define the maximum number of opacity species your program will handle
#define MAX_OPACITY_SPECIES 100

// Define the maximum length for species names
#define MAX_NAME_LENGTH 50


/* --- Opacity structure -------------------------------------------- */
struct Opac {
    int NP; // Number of pressure points
    int NT; // Number of temperature points
    double *T; // Array for temperature values
    double *P; // Array for pressure values
    double *Plog10; // Array for log10 of pressure values
    double ***kappa; // 3D array for opacity values [species][pressure][temperature]
    double **abundance; // 2D array for abundance [species][pressure or temperature]
    char **speciesNames; // Array of pointers to species name strings
    int numSpecies; // Number of opacity species
};

/* --- Chemistry structure ------------------------------------------ */
struct Chem {
    double ***speciesData; // 3D array of species data [species][pressure][temperature]
    double *T;  // Array for temperature values
    double *P;  // Array for pressure values
    char **species; // Array of species names
    int numSpecies; // Number of chemical species
};

// Structure to hold species names
typedef struct {
  char names[MAX_OPACITY_SPECIES][MAX_NAME_LENGTH];
  int count;
} SpeciesList;

/* --- Global variables --------------------------------------------- */
// If you are using these structures as global variables, declare them here.
// Otherwise, manage them locally within your functions as needed.
extern struct Opac *opac;
extern struct Chem chem;

/* Function prototypes */
void TotalOpac(); // No argument needed as we're using global variables
void FreeOpacTable(struct Opac* opac);
void ReadChemTable(void);
void FreeChemTable(void);

#endif /* !__OPAC_H__ */
