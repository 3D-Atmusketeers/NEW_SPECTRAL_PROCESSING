/* ------- file: -------------------------- opac.h ------------------ */

#ifndef __OPAC_H__
#define __OPAC_H__


#define MAX_OPACITY_SPECIES 30 // Maximum number of species names
#define MAX_NAME_LENGTH 200

/* Variables for opacity
 *
 * Sara Seager
 * Last modified: Eliza Miller-Ricci June 25, 2007
 */

/* --- Opacity structure -------------------------------------------- */
struct Opac {
    // Existing fields...
    int NP;
    int NT;
    double *T;
    double *P;
    double *Plog10;
    double ***kappa;
    double **abundance;
    char speciesName[MAX_NAME_LENGTH];
};



/* --- Chemistry structure ------------------------------------------ */

// Define the structure for chemical data
struct Chem {
    double ***speciesData; // Array of species data
    double *T;  // Temperature array
    double *P;  // Pressure array
    char **species; // Array of species names
    int numSpecies; // Number of species
};


// Declare the Chem structure
extern struct Chem chem;

/* --- Associated function prototypes ------------------------------- */

void ReadOpacTable(struct Opac *opac, const char *filename, const char *speciesName);
void FreeOpacTable(struct Opac* opac);
void TotalOpac();
void ReadChemTable();
void FreeTP();

#endif /* !__OPAC_H__ */

/* ------- end ---------------------------- opac.h ------------------ */
