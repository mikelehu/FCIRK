/*-----------------------------------------------------------------------------*/
/*									       */
/*                       Problems.h	  		       		       */
/*									       */
/* ----------------------------------------------------------------------------*/


#include <def.h>
#include <Kepler.h>


/*    Functions to recovery the barycentric coordinates from                       */
/*    canonical heliocentric coordinates                                           */

void ChangeHeltoBar_EMB (int neqH, int nbodyH, val_type *u, solution *U,
                         val_type *Gm, val_type *mu);
                         
void ChangeHeltoBar_Moon (int neqH, int nbodyH, val_type *u, solution *U, 
                          val_type *Gm, val_type *mu);                         


/******************** Hamiltonian ****************/

__float128 HamNbody (int neq,solution *u,parameters *params);
__float128 Ham_K (int nbody,solution *u, parameters *params);


/******************** 1-Problem: Earth+Moon one body****************/

__float128 RR1 (int neq,solution *u,parameters *params);
void DRR1 (int neq, val_type t,val_type *u,val_type *dR,parameters *params);
void DRR1_high (int neq, highprec t,highprec *u,highprec *dR,parameters_high *params);
void KComp1 (int nbody, ode_sys *system, Pkepler_sys *Pkepler);

/******************** Perturbations RR2: Moon as separate body************/

__float128 RR2 (int neq,solution *u,parameters *params);
void DRR2 (int neq, val_type t,val_type *u,val_type *dR,parameters *params);
void DRR2_high (int neq, highprec t,highprec *u,highprec *dR,parameters_high *params);
void KComp2 (int nbody, ode_sys *system, Pkepler_sys *Pkepler);
void Mysubstract (val_type *q, val_type *Dq, val_type *subs);
void Mysubstract_high (highprec *q, highprec *Dq, highprec *sub);

/******************integration based on Kepler flow's change variables******/

void NbodyGFcn (int neq, val_type t,val_type *u,val_type *dR,parameters *params);
void NbodyGFcn_high (int neq, highprec t,highprec *u,highprec *dR,parameters_high *params);
void NbodyOde (int neq, val_type t, val_type ttau,val_type *u,val_type *f
              ,parameters *params,tcache_vars *cache_vars);
void NbodyOde_high (int neq, highprec t,highprec ttau,highprec *u,highprec *f,
                    parameters_high *params,tcache_vars_high *cache_vars_high);


