#include <stdio.h>
#include <stdlib.h>
#include "input.h"
#include "opac.h"
#include "atmos.h"

/* --- Global variables ------------------------------------------ */
struct Atmos atmos;
struct Opac* opac;
extern struct Chem chem; // Assuming chem is declared in opac.h and defined in totalopac.c or another related source file.

/* --- Function prototypes --------------------------------------- */
double PHASE;

void TotalOpac(); // Updated to reflect no arguments as per totalopac.c implementation.
void ReadTP_3D();
void RT_Emit_3D(double PHASE, struct Opac* opac);
void FreeChemTable(void); // Assuming declaration in opac.h.

/* ------- begin ---------------- main --------------------------- */
int main() {
    PHASE = 0.0;
    double phase_step = 360. / N_PHASE;

    // Initialize chemical and opacity data structures
    ReadChemTable(); // This prepares the chem structure with species data.

    int i;
    for (i = 0; i < N_PHASE; i++) {
        printf("Phase: %06.2f\n", PHASE);

        // Calculate total opacity based on current chemical data and other factors.
        TotalOpac(); // Updated to no longer require passing an Opac structure, as opac handling is internal to the function.

        printf("TotalOpac done\n");

        ReadTP_3D();
        printf("ReadTP done\n");

        // Assuming RT_Emit_3D requires an Opac structure, it must be updated to handle the global opac or an array of Opac structures.
        RT_Emit_3D(PHASE, opac); // Needs adjustment based on how opac is managed after integrating totalopac.c code.
        printf("RT_Emit done\n");

        PHASE += phase_step;
        printf("\n");
    }

    // Assuming that FreeChemTable is sufficient for cleanup as TotalOpac and related functions manage Opac structure(s).
    FreeChemTable(); // Cleanup for the Chem structure.

    // Additional cleanup might be required here for the Opac structures if they are allocated dynamically.

    return 0;
}
