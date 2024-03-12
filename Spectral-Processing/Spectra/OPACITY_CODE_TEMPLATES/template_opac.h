/* ------- file: -------------------------- opac.h ------------------ */

#ifndef __OPAC_H__
#define __OPAC_H__

/* --- Opacity structure -------------------------------------------- */
struct Opac { 
  int NP, NT;
  double ***kappa, ***alpha, **kappa_pl, **kappa_ross, **kappa_mean, 
    **abundance, **mu;
  double *P, *Plog10, *T;
};

/* --- Chemistry structure ------------------------------------------ */

struct Chem {
  double  **total,
          **CH4,
          **PH3,
          **K,
          **CrH,
          **He,
          **H2S,
          **Ca+,
          **VO,
          **SiH4,
          **N,
          **Fe+,
          **H2,
          **TiO,
          **Fe,
          **el,
          **CO2,
          **CaH,
          **SiH,
          **SiO2,
          **C,
          **H2O,
          **SH,
          **SiO,
          **CaO,
          **AlH1,
          **Na,
          **O3,
          **AlO,
          **CO,
          **C2H2,
          **N2,
          **O2,
          **NaH,
          **Ca,
          **NH3,
          **MgH,
          **FeH,
          **NO,
          **SO2,
          **HF,
          **O,
          **HO,
          **Mg,
          **H,
          **HCN,
          **map_abund,
          **rho_k_T,
          **mu;
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
