/* ------- file: -------------------------- opac.h ------------------ */

#ifndef __OPAC_H__
#define __OPAC_H__

/* Variables for opacity
 *
 * Sara Seager
 * Last modified: Eliza Miller-Ricci June 25, 2007
 */

/* --- Opacity structure -------------------------------------------- */
struct Opac { 
  int NP, NT;
  double ***kappa, ***alpha, **kappa_pl, **kappa_ross, **kappa_mean, 
    **abundance, **mu;
  double *P, *Plog10, *T;
};

/* --- Chemistry structure ------------------------------------------ */

struct Chem {
  double **total,
                  **E,
                  **H,
                  **H2,
                  **He,
                  **C2H2,
                  **CH4,
                  **CO,
                  **CO2,
                  **H2O,
                  **H2S,
                  **HCN,
                  **K,
                  **Na,
                  **NH3,
                  **PH3,
                  **TiO,
                  **VO,
                  **FeH,
                  **Ar,
                  **N2,
                  **O2,
                  **map_abund,
                  **rho_k_T, **mu;
  double *P, *T;
};

/* --- Associated function prototypes ------------------------------- */

void ReadOpacTable(struct Opac opac, char *filename);
void FreeOpacTable(struct Opac opac);
void TotalOpac();
void ReadChemTable();
void FreeTP();

#endif /* !__OPAC_H__ */

/* ------- end ---------------------------- opac.h ------------------ */
