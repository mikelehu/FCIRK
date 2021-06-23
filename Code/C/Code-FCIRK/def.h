/*----------------------------------------------------------------------------*/
/*									      */
/*                                def.h                                       */
/*					                                      */
/* ---------------------------------------------------------------------------*/

#include <stdbool.h>
#include <float.h>

/* ---------------------------------------------------------------------------*/
/*									      */
/*	Parameters							      */
/*									      */
/* ---------------------------------------------------------------------------*/

#ifndef DEFH_
#define DEFH_

#define INF DBL_MAX
#define MAXIT 100		// Maximum number of fixed point iterations

#define SUCCESS 0
#define FAIL -1
#define STRMAX 256		// Filename maximum string
#define RTOL pow(10.,-12);	// pow(2.,-40);
#define ATOL pow(10.,-12);
#define MAXPARAM 30		// Number of maximum parameters
#define NU 1000                 // Threshold for adaptive step 100
#define KMAX 10                 // Max number of substeps (adaptive) 

#define DIR_COEFF "../../../Code/C/CoefficientsData/" // Path for Coefficients 


/* ---------------------------------------------------------------------------*/
/*									      */
/*	PREC: 								      */
/*	  1=Double - Long Double					      */
/*	  2=Double - quadruple						      */
/*	  3=Long Double - quadruple					      */
/*	  11=quadruple							      */
/*									      */
/* ---------------------------------------------------------------------------*/


#if PREC ==1 // double-long double
//
typedef double val_type;
typedef long double highprec;
#define TOLKEPLER POW(10,-20)
#define TOLKEPLER_high POW(10,-30)
#define POW(x, y)      pow(x, y)
#define SQRT(x)        sqrt(x)
#define EXP(x)         exp(x)
#define LOG(x)         log(x)
#define SIN(x)         sin(x)
#define COS(x)         cos(x)
#define SINCOS(x,s,c)  sincos(x,s,c)
#define TAN(x)         tan(x)
#define ATAN2(x)       atan2(x)
#define FABS(x)	fabs(x)
#define FMAX(x,y)      fmax(x,y)
//
#define POW_high(x, y)      powl(x, y)
#define SQRT_high(x)        sqrtl(x)
#define EXP_high(x)         expl(x)
#define LOG_high(x)         logl(x)
#define SIN_high(x)         sinl(x)
#define SINCOS_high(x,s,c)  sincosl(x,s,c)
#define COS_high(x)         cosl(x)
#define TAN_high(x)         tanl(x)
#define ATAN2_high(x)       atan2l(x)
#define FABS_high(x)        fabsl(x)
#define FMAX_high(x,y)      fmaxl(x,y)
//
#elif PREC == 2 // double-quadruple
//
typedef double val_type;
typedef __float128 highprec;
#define TOLKEPLER pow(10,-20)
#define TOLKEPLER_high pow(10,-40)
#define POW(x, y)      pow(x, y)
#define SQRT(x)        sqrt(x)
#define EXP(x)         exp(x)
#define LOG(x)         log(x)
#define SIN(x)         sin(x)
#define SINCOS(x,s,c)  sincos(x,s,c)
#define COS(x)         cos(x)
#define TAN(x)         tan(x)
#define ATAN2(x)       atan2(x)
#define FABS(x)	       fabs(x)
#define FMAX(x,y)      fmax(x,y)
//
#define POW_high(x, y)      powq(x, y)
#define SQRT_high(x)        sqrtq(x)
#define EXP_high(x)         expq(x)
#define LOG_high(x)         logq(x)
#define SIN_high(x)         sinq(x)
#define COS_high(x)         cosq(x)
#define SINCOS_high(x,s,c)  sincosq(x,s,c)
#define TAN_high(x)         tanq(x)
#define ATAN2_high(x)       atan2q(x)
#define FABS_high(x)        fabsq(x)
#define FMAX_high(x,y)      fmaxq(x,y)
//
#elif PREC ==3  // long double-quadruple
//
typedef long double val_type;
typedef __float128 highprec;
#define TOLKEPLER pow(10,-25)
#define TOLKEPLER_high pow(10,-25)    
#define POW(x, y)      powl(x, y)
#define SQRT(x)        sqrtl(x)
#define EXP(x)         expl(x)
#define LOG(x)         logl(x)
#define SIN(x)         sinl(x)
#define COS(x)         cosl(x)
#define SINCOS(x,s,c)  sincosl(x,s,c)
#define TAN(x)         tanl(x)
#define ATAN2(x)       atan2l(x)
#define FABS(x)	       fabsl(x)
#define FMAX(x,y)      fmaxl(x,y)
//
#define POW_high(x, y)      powq(x, y)
#define SQRT_high(x)        sqrtq(x)
#define EXP_high(x)         expq(x)
#define LOG_high(x)         logq(x)
#define SIN_high(x)         sinq(x)
#define COS_high(x)         cosq(x)
#define SINCOS_high(x,s,c)  sincosq(x,s,c)
#define TAN_high(x)         tanq(x)
#define ATAN2_high(x)       atan2q(x)
#define FABS_high(x)        fabsq(x)
#define FMAX_high(x,y)      fmaxq(x,y)
//
#elif PREC == 11           // quadprecision 
//
typedef __float128 val_type;
typedef __float128 highprec;
#define TOLKEPLER pow(10,-40)
#define TOLKEPLER_high pow(10,-40)
#define POW(x, y)      powq(x, y)
#define SQRT(x)        sqrtq(x)
#define EXP(x)         expq(x)
#define LOG(x)         logq(x)
#define SIN(x)         sinq(x)
#define COS(x)         cosq(x)
#define SINCOS(x,s,c)  sincosq(x,s,c)
#define TAN(x)         tanq(x)
#define ATAN2(x)       atan2q(x)
#define FABS(x)	       fabsq(x)
#define FMAX(x,y)      fmaxq(x,y)
//
#define POW_high(x, y)      powq(x, y)
#define SQRT_high(x)        sqrtq(x)
#define EXP_high(x)         expq(x)
#define LOG_high(x)         logq(x)
#define SIN_high(x)         sinq(x)
#define COS_high(x)         cosq(x)
#define SINCOS_high(x,s,c)  sincosq(x,s,c)
#define TAN_high(x)         tanq(x)
#define ATAN2_high(x)       atan2q(x)
#define FABS_high(x)        fabsq(x)
#define FMAX_high(x,y)      fmaxq(x,y)
#endif


/* ---------------------------------------------------------------------------*/
/*									      */
/*	General definitions					       	      */
/*								              */
/* ---------------------------------------------------------------------------*/


typedef struct tcoeffs
   {
     int ns;
     val_type *c,*b,*a;	 	 // c,b,a coefficients.
     val_type *m;	 	 // mij=aij/bj and mij+mji-1=0.
     val_type *hc;       	 // hc=h*c.
     val_type *hb;       	 // hb=h*b.
     val_type *nu;               // interpolate coefficients (a*/bj).
     val_type *ttau;             // ttau=hc-h/2
     val_type *d;		 // coefficients for adaptive step

   } tcoeffs;


typedef struct tcoeffs_h
   {
     int ns;
     highprec *c,*b,*a;	 	 // c,b,a coefficients.
     highprec *m;	 	 // mij=aij/bj and mij+mji-1=0.
     highprec *hc;       	 // hc=h*c.
     highprec *hb;       	 // hb=h*b.
     highprec *nu;               // interpolate coefficients (a*/bj).
     highprec *ttau;             // ttau=hc-h/2
     highprec *d;		 // coefficients for adaptive step

   } tcoeffs_h;


typedef struct tmethod 
  {
    tcoeffs       coeffs;
    tcoeffs_h     coeffs_h;

  } tmethod;

typedef struct solution
  {
     highprec *uu;
     val_type *uul;

  } solution;


typedef struct toptions
  {
     val_type *rtol,*atol;
     int sampling;
     bool adaptive;
     char filename[STRMAX];     // Output filename.
     void (*TheOutput)();       // Output function.

   } toptions;


typedef struct Pkepler_sys
   {
     __float128 (*RR)();     // hamiltonian perturbation part.
     void (*DRR)();          // {dRR/dq, dRR/dp}
     void (*DK)();           // {Keplerian part}
     void (*KComputation)(); 
     int keplerkop;
     val_type *K;           // OdeKepler real parameters
     val_type *K2;          
     __float128 *Khigh;
     val_type *Mu;
     __float128 *Muhigh;               

   }  Pkepler_sys;

typedef struct Pkepler_sys_high
   {
     __float128 (*RR)();     // hamiltonian perturbation part.
     void (*DRR)();          // {dRR/dq, dRR/dp}
     void (*DK)();           // {Keplerian part} 
     int keplerkop;
     highprec *K;           // OdeKepler real parameters
     highprec *K2;          
     highprec *Mu;             
   }  Pkepler_sys_high;

typedef struct parameters
   {
     val_type *rpar;	 // Variables for specifying odefun real parameters.
     __float128 *rparhigh; 
     int numrpar;	 // Number of real parameters.
     int *ipar;          // Variables for specifying odefun integer parameters.
     int numipar;	 // Number of int parametes.
     Pkepler_sys Pkepler;

   } parameters;

typedef struct parameters_high
   {
     highprec *rpar;	 // Variables for specifying odefun real parameters.
     int numrpar;	 // Number of real parameters.
     int *ipar;          // Variables for specifying odefun integer parameters.
     int numipar;	 // Number of int parametes.
     Pkepler_sys_high Pkepler;

   } parameters_high;


typedef struct ode_sys
   {
     int neq;		        // number of equations.
     void (*f)();	        // odefun.
     __float128 (*ham)();	// hamiltonian
     parameters params;
     void (*ProjFun)();      
     void (*StartFun)();     
     int ode;       

    } ode_sys;


typedef struct ode_sys_high
   {
     int neq;		        // number of equations.
     void (*f)();	        // odefun.
     __float128 (*ham)();	// hamiltonian
     parameters_high params;
     int ode;       

    } ode_sys_high;


typedef struct tode_sys
   {
     ode_sys system;
     ode_sys_high system_h;

   } tode_sys;

typedef struct tcache_stat
    {
    int execution;            // SUCCESS or FAIL.
    int convergence;	      // SUCCESS or FAIL.
    bool laststep;     

    __float128 E0;     	      // Initial Energy
    double MaxDE; 	      // MaxDE=Abs(Ei-E0/E0)
    bool interpolate;        //  interpolate=true or interpolate=false

    /* statistics*/
    int stepcount;
    int hstepcount;
    int itcount;
    long int totitcount;
    long int totitcountzero;	
    int maxitcount;
    int itzero;
    int fcn;
    int nout;                // number of output values.

    } tcache_stat;


typedef struct tcache_vars
    {

    val_type *z,*li,*fz;
    val_type *zold;
    val_type *DMin;
    val_type *delta, *avg_delta;          

    } tcache_vars;

typedef struct tcache_vars_high
 {
     highprec *z,*li,*fz;
     highprec *zold; 
     highprec *DMin;  

 } tcache_vars_high;


typedef struct tcache
 {
     tcache_stat cache_stat;
     tcache_vars cache_vars;
     tcache_vars_high cache_vars_h;

 } tcache;

typedef struct flowaux
   {
    val_type r0;
    val_type eta;
    val_type alpha;
    val_type beta;
    val_type gamma;
    val_type zeta;
    val_type X;
    val_type GG[3];
    val_type r;
    val_type rinv;
    val_type b[4];

   } flowaux;

typedef struct flowaux_high
   {
    highprec r0;
    highprec eta;
    highprec alpha;
    highprec beta;
    highprec gamma;
    highprec zeta;
    highprec X;
    highprec GG[3];
    highprec r;
    highprec rinv;
    highprec b[4];

   } flowaux_high;

#endif /*DEFH_*/
