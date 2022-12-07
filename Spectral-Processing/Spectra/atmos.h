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
    ***incident_frac,
    ***aero_tau_pre_qext_1, ***sw_asym_1, ***sw_pi0_1,
    ***aero_tau_pre_qext_2, ***sw_asym_2, ***sw_pi0_2,
    ***aero_tau_pre_qext_3, ***sw_asym_3, ***sw_pi0_3,
    ***aero_tau_pre_qext_4, ***sw_asym_4, ***sw_pi0_4,
    ***aero_tau_pre_qext_5, ***sw_asym_5, ***sw_pi0_5,
    ***aero_tau_pre_qext_6, ***sw_asym_6, ***sw_pi0_6,
    ***aero_tau_pre_qext_7, ***sw_asym_7, ***sw_pi0_7,
    ***aero_tau_pre_qext_8, ***sw_asym_8, ***sw_pi0_8,
    ***aero_tau_pre_qext_9, ***sw_asym_9, ***sw_pi0_9,
    ***aero_tau_pre_qext_10, ***sw_asym_10, ***sw_pi0_10,
    ***aero_tau_pre_qext_11, ***sw_asym_11, ***sw_pi0_11,
    ***aero_tau_pre_qext_12, ***sw_asym_12, ***sw_pi0_12,
    ***aero_tau_pre_qext_13, ***sw_asym_13, ***sw_pi0_13,
    ***aero_tau_haze, ***tau_asym, ***tau_pi0;
  double startlambda, endlambda;




};

/* --- Associated function prototypes ------------------------------- */

void ReadTP();
void FreeTP();

#endif /* !__ATMOS_H__ */

/* ------- end ---------------------------- atmos.h ----------------- */
