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
( ode_sys *system, gauss_method *gsmethod, solver_stat *thestatptr
);

val_type NormalizedDistance
( int neq, int ns, toptions *options,
  val_type *z, val_type *zold
);


void Fixed_point_Step
( ode_sys *system,   solution *u,
  val_type tn,  val_type h,
  toptions *options,  gauss_method *method,
 solver_stat *thestatptr
);

int General_FP_It			
( ode_sys *system,  solution *u,  val_type tn,
  val_type h,  gauss_method *method,
  solver_stat *thestatptr,
  int *D0, bool *iter0, val_type *DMin
);


void Summation
( gauss_method *gsmethod,
 solution *u,
  ode_sys *system,  toptions *options,
  solver_stat *thestatptr
);


void Main_FCIRK            
(val_type t0, val_type t1, val_type h,
  gauss_method *gsmethod, solution *u,
  ode_sys *system, toptions *options,
 solver_stat *thestatptr
);


