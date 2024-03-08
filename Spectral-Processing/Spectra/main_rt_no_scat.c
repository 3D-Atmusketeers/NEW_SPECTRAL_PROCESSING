#include <stdio.h>
#include "input.h"
#include "opac.h"
#include "atmos.h"
#include <stdlib.h>

/* --- Global variables ------------------------------------------ */

struct Atmos atmos;
struct Opac* opac;

/* --- Function prototypes --------------------------------------- */

double PHASE;

void TotalOpac(struct Opac* opac); 
void ReadTP_3D();
void RT_Emit_3D(double PHASE, struct Opac* opac);

/* ------- begin ---------------- main --------------------------- */

int main() {
    PHASE = 0.0;
    double phase_step = 360. / N_PHASE;
    
    // Dynamically allocate the Opac structure.
    opac = (struct Opac*)malloc(sizeof(struct Opac));
    if (opac == NULL) {
        fprintf(stderr, "Failed to allocate memory for Opac structure.\n");
        return EXIT_FAILURE;
    }
    
    int i;
    
    for(i = 0; i < N_PHASE; i++) {
        printf("Phase: %06.2f\n", PHASE);
        
        // Pass the dynamically allocated Opac structure.
        TotalOpac(opac); 

        printf("TotalOpac done\n");

        // Assume this does not need updating.
        ReadTP_3D(); 

        printf("ReadTP done\n");

        // Pass the Opac structure as needed.
        RT_Emit_3D(PHASE, opac); 
        printf("RT_Emit done\n");
        
        PHASE += phase_step;
        printf("\n");
    }

    // Free the dynamically allocated Opac structure after use.
    FreeOpacTable(opac); // Assuming you have implemented this function.
    free(opac);
    opac = NULL;

    return 0;
}
