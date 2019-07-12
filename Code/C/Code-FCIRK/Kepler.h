/*------------------------------------------------------------------------------*/
/*										*/
/*                                Kepler.h					*/
/*										*/
/* -----------------------------------------------------------------------------*/

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h>
#include <math.h>
#include <def.h>

int sign
 (val_type a
);

int sign_high 
 (highprec a
);

void KeplerFlowAll
(int neq,int keplerkop, solution *u, 
 val_type h, parameters *params
);


void KeplerFlowAll2
(int neq,int keplerkop, val_type *u,
 val_type h, parameters *params
);


void KeplerFlowAll_high
(int neq,int keplerkop, solution *u, 
 highprec h, parameters *params
);

void KeplerFlow
(val_type mu, val_type *q,val_type *v,
 val_type *dq,val_type *dv,val_type dt
);

void KeplerFlow_high
(highprec mu, highprec *q,highprec *v,
 highprec *dq,highprec *dv,highprec dt
);

void KeplerFlowMIXED_high 
(highprec mu, highprec *q,highprec *v,
 highprec *dq,highprec *dv,highprec dt
);

void KeplerSolveE
(val_type r0, val_type gamma, val_type eta, 
 val_type beta, val_type k, val_type dt,
 val_type *X, val_type *G
); 


void KeplerSolveE_high
(highprec r0, highprec gamma, highprec eta, 
 highprec beta, highprec k, highprec dt,
 highprec *X, highprec *G
); 


void KeplerSolveITER_high
(highprec r0, highprec gamma, highprec eta, 
 highprec beta, highprec k, highprec dt,
 highprec *X, highprec *G
);

void  KeplerFlowGFcn 
(void GFcn(), int neq, val_type t,
 val_type *U, int keplerkop,
 val_type *mu, parameters *params, val_type dt,
 val_type *G
);


void KeplerFlowGen 
(val_type dt, val_type *Q, val_type *V, val_type mu,
 val_type *qq, val_type *vv, flowaux *aux
);


void KeplerFlowGFcnAux
(val_type dt,val_type *Q, val_type *V, 
 val_type *gq,val_type *gv,
 flowaux *aux, val_type mu, val_type *Gq, val_type *Gv
);


void ProjFun 
(int neq, 
 val_type t, val_type h, solution *w, parameters *params
);

void ProjFun_high
(int neq, val_type t, val_type h, solution *w, parameters *params
);

void StartFun
(int neq, val_type t, val_type h, solution *w, parameters *params
);

void StartFun_high
(int neq, val_type t, val_type h, solution *w, parameters *params
);

void OutputFun 
( ode_sys *system,  gauss_method *method,
  val_type t, val_type h, solution *w,
  solver_stat *thestatptr,
  parameters *params, toptions *options,FILE *loga
);

void OutputFun_high 
( ode_sys *system,  gauss_method *method,
  val_type t, val_type h, solution *w,
  solver_stat *thestatptr,
  parameters *params, toptions *options,FILE *loga
);


