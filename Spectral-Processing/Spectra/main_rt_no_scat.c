
#include <stdio.h>
#include "input.h"
#include "opac.h"
#include "atmos.h"

/* --- Global variables ------------------------------------------ */

struct Atmos atmos;
struct Opac opac;

/* --- Function prototypes --------------------------------------- */

double PHASE;

void TotalOpac();
void ReadTP_3D();
void RT_Emit_3D(double PHASE);

/* ------- begin ---------------- main --------------------------- */


int main(){
    
    PHASE = 0.0;
    double phase_step = 360. / N_PHASE;
    
    int i;
    
    for(i=0; i<N_PHASE; i++){
        printf("Phase: %06.2f\n", PHASE);
        
        TotalOpac();

        printf("TotalOpac done\n");

        ReadTP_3D();

        printf("ReadTP done\n");

        RT_Emit_3D(PHASE);
        printf("RT_Emit done\n");
        
        PHASE += phase_step;
        printf("\n");
        
    }
  return 0;
}

