/*-----------------------------------------------------------------------------*/
/*									       */
/*                       Problems.h	  		       		       */
/*									       */
/* ----------------------------------------------------------------------------*/


#include <def.h>
#include <Kepler.h>


/******************** Keplerian part*****************************************/
__float128 HamKepler (int nbody,solution *u, parameters *params);
void OdeKepler (int neq, val_type t,val_type *u,val_type *f,parameters *params);

/******************* Perturbation based odefun *****************************/

/* Generic Odefun: Dkepler+RR (IRKFP, CO1035)*/
void Ode3 (int neq, val_type t,val_type *u,val_type *f,parameters *params);
__float128 Ham3 (int neq,solution *u,parameters *params);

/******************integration based on Kepler flow's change variables******/

void NbodyGFcn (int neq, val_type t,val_type *u,val_type *dR,parameters *params);
void NbodyOde (int neq, val_type t, val_type ttau,val_type *u,val_type *f,parameters *params);

/******************** Perturbation RR1: Earth+Moon one body****************/
/******* Equations corresponding to Canonical Heliocentric Coordinates ****/

/* EDO Perturbations */
__float128 RR1 (int neq,solution *u,parameters *params);

/* Hamiltonian part*/
void DRR1 (int neq, val_type t,val_type *u,val_type *dR,parameters *params);

/* Computation of k values*/
void KComp1 (int nbody, ode_sys *system, Pkepler_sys *Pkepler);

/******************** Perturbations RR2: Moon as separate body************/

/* EDO Perturbations */
__float128 RR2 (int neq,solution *u,parameters *params);
void mysubstract (val_type *q, val_type *Dq, val_type *subs);

/* Hamiltonian part*/
void DRR2 (int neq, val_type t,val_type *u,val_type *dR,parameters *params);

/* Computation of k values*/
void KComp2 (int nbody, ode_sys *system, Pkepler_sys *Pkepler);



