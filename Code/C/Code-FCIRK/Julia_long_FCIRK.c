/*----------------------------------------------------------------------------*/
/*									      */
/*                         Julia_long_FCIRK.c	       			      */
/*									      */
/* ---------------------------------------------------------------------------*/

#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <def.h>
#include <stdio.h>
#include <Problems.h>
#include <Common-FCIRK.h>
#include <GaussCoefficients.h>
#include <Kepler.h>
#include <time.h>
#include <stdbool.h>
#include <quadmath.h>
#include <gmp.h>
#include <mpfr.h>



/* Global variables */

int thread_count;


int longFCIRK (int Nkepler,int Moreq, 
               int ns, double t0, double t1,
               mpfr_ptr *u0, int ulen,
               double h,
               mpfr_ptr *k, int klen,
               mpfr_ptr *rpar,int rlen,
               int *ipar,int ilen,
               const char *myfilename,
               int sampling,
               int ode,
               int threads,
               bool adaptive,
               bool mixed,
               bool monitoring_err,
               int nrmbits,
               double *result_array)

{

#include <Julia_Source_FCIRK.h>

}



