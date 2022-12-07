/* ------- file: -------------------------- atmos.h ----------------- */


#ifndef __ATMOS_H__
#define __ATMOS_H__

/* Defines Atmos structure.
 *
 * Eliza Miller-Ricci
 * Last modified: 
 */

/* --- Structure defines atmosphere --------------------------------- */
 
struct Atmos {

  int Nlambda;
  double **kappa_nu, **tau_nu, **dtau_nu;
  double *lambda, *nu, *tau, *T, *P, *J, *H, *kappa_m, *kappa_kg, *rho, *ds, 
    *dtau, *mu, *lat, *lon, *alt, ***T_3d, ***P_3d, ***vel_ew, ***vel_ns, ***vel_ve,
    ***aero_sw_tau_1, ***sw_asym_1, ***sw_pi0_1,
    ***aero_sw_tau_2, ***sw_asym_2, ***sw_pi0_2,
    ***aero_sw_tau_3, ***sw_asym_3, ***sw_pi0_3,
    ***aero_sw_tau_4, ***sw_asym_4, ***sw_pi0_4;
  double startlambda, endlambda;

};

/* --- Associated function prototypes ------------------------------- */

void ReadTP();
void FreeTP();

#endif /* !__ATMOS_H__ */

/* ------- end ---------------------------- atmos.h ----------------- */
