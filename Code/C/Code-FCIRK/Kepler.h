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


void KeplerSolveE
(val_type r0, val_type eta, val_type zeta,  
 val_type beta, val_type k, val_type t,
 val_type *X, val_type *G, int *iter
);

void KeplerSolveE_high
(highprec r0, highprec eta, highprec zeta, 
 highprec beta, highprec k, highprec t,
 highprec *X, highprec *G, int *iter
 ); 


void  KeplerFlowGFcn 
(void GFcn(), int neq, val_type t,
 val_type *U, int keplerkop,
 val_type *k, parameters *params, val_type dt,
 val_type *G
);


void  KeplerFlowGFcn_high 
(void GFcn(), int neq, highprec t,
 highprec *U, int keplerkop,
 highprec *k, parameters_high *params, highprec dt,
 highprec *G
);

void KeplerFlowGen 
(val_type dt, val_type *Q, val_type *V, val_type mu,
 val_type *qq, val_type *vv, flowaux *aux
);


void KeplerFlowGen_high
(highprec dt, highprec *Q, highprec *V, highprec k,
 highprec *qq, highprec *vv, flowaux_high *aux
);

void KeplerFlowGFcnAux
(val_type dt,val_type *Q, val_type *V, 
 val_type *gq,val_type *gv,
 flowaux *aux, val_type mu, val_type *Gq, val_type *Gv
);

void KeplerFlowGFcnAux_high
(highprec dt,highprec *Q, highprec *V, 
 highprec *gq, highprec *gv,
 flowaux_high *aux, highprec k, highprec *Gq, highprec *Gv
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
( ode_sys *system,  tcoeffs *method,
  val_type t, val_type h, solution *w,
  tcache_stat *cache_stat,
  parameters *params, toptions *options,FILE *loga
);

void OutputFun_high 
( ode_sys *system,  tcoeffs *method,
  val_type t, val_type h, solution *w,
  tcache_stat *cache_stat,
  parameters *params, toptions *options,FILE *loga
);


