/*-----------------------------------------------------------------------------*/
/*									       */
/*                        Common-FCIRK.h	         		       */
/*									       */
/* ----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <def.h>
#include <sys/stat.h>
#include <Problems.h>
#include <GaussCoefficients.h>
#include <Kepler.h>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <quadmath.h>



void InitStat
(ode_sys *system, tcoeffs *gsmethod,
 tcache_stat *cache_stat, tcache_vars *cache_vars
);


void add2 
(val_type x, val_type xx, val_type y, val_type yy,
           val_type *z, val_type *zz
);

void InitStat_high
(ode_sys_high *system, tcoeffs_h *gsmethod,
 tcache_vars_high *cache_vars
);


val_type Rmdigits
(val_type x,val_type r
);

highprec Rmdigits_high 
(highprec x,highprec r
);

val_type NormalizedDistance
( int neq, int ns, toptions *options,
  val_type *z, val_type *zold
);

val_type NormalizedDistance_high
( int neq, int ns,
  toptions *options,  highprec *z,
  highprec *zold
);

void IRKstep_fixed
( ode_sys *system,   solution *u,
  val_type tn,  int ii,val_type h,
  toptions *options,  tcoeffs *method,
  tcache_stat *cache_stat, tcache_vars *cache_vars
);

void IRKstep_fixed_high
( ode_sys_high *system,  solution *u,
  highprec tn, int ii, highprec h,
  toptions *options,  tcoeffs_h *method,
  tcache_stat *cache_stat, tcache_vars_high *cache_vars
);


void IRKstep_adaptive
(tode_sys *ode_system, solution *u,
 val_type tn, int ii, val_type h,
 toptions *options,  tmethod *method, 
 tcache_stat *cache_stat, 
 tcache_vars *cache_vars, tcache_vars_high *cache_vars_high
);


int FP_Iteration			
( ode_sys *system,  solution *u,  val_type tn,
  int ii, val_type h, 
  toptions *options, tcoeffs *method,
  tcache_stat *cache_stat, tcache_vars *cache_vars,
  int *D0, bool *iter0
);

int FP_Iteration_high
( ode_sys_high *system,  solution *u,  highprec tn,
  int ii, highprec h, 
  toptions *options, tcoeffs_h *method,
  tcache_stat *cache_stat, tcache_vars_high *cache_vars,
  int *D0, bool *iter0
);

void Summation
( tcoeffs *gsmethod,
  solution *u,
  ode_sys *system,  toptions *options,
  tcache_vars *cache_vars
);


void Summation_high
( tcoeffs_h *gsmethod,
  solution *u, ode_sys_high *system,
  toptions *options, tcache_vars_high *cache_vars
);


void deltafun
( tcoeffs *gsmethod, ode_sys *system,
  tcache_vars *cache_vars, val_type *delta
);

bool Ordinary_stepQ 
( ode_sys *system,
  val_type h,
  int k,
  tcache_stat *cache_stat, 
  tcache_vars *cache_vars
);

int Num_steps 
( tcoeffs *gsmethod,
  ode_sys *system,
  tcache_stat *cache_stat,
  tcache_vars *cache_vars
);


void Main_FCIRK
( val_type t0, val_type t1, val_type h, 
  tmethod *method, solution *u, 
  tode_sys *ode_system, toptions *options,
  tcache *cache
);









