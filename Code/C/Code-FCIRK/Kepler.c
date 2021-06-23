/*------------------------------------------------------------------------------*/
/*										*/
/*     Kepler.c									*/
/*										*/
/*	Functions:								*/
/*										*/
/*	   sign									*/
/*	   sign_high								*/
/*										*/
/*         KeplerFlowAll()							*/
/*         KeplerFlowAll2()							*/
/*	   KeplerFlowAll_high							*/
/*	   KeplerFlow								*/
/*	   KeplerFlow_high							*/
/*	   KeplerSolveE								*/
/*	   KeplerSolveE_high							*/
/*	   KeplerFlowGFcn							*/
/*         KeplerFlowGen							*/
/*         KeplerFlowGFcnAux							*/
/*	   ProjFun								*/
/*	   ProjFun_high								*/
/*	   StartFun								*/
/*	   StartFun_high						        */
/*	   OutputFun								*/
/*	   OutputFun_high							*/
/*										*/
/*	Related code files (source):						*/
/*										*/
/*	   source_KeplerFlowAll.c						*/
/*	   source_KeplerFlow.c							*/
/*	   source_KeplerSolveE.c						*/
/*	   source_ProjFun.c							*/
/*	   source_StartFun.c							*/
/*	   source_OutputFun.c							*/
/*										*/
/* -----------------------------------------------------------------------------*/


#include <Kepler.h>
#include <quadmath.h>
#include <Common-FCIRK.h>

void sincos(double x, double *sin, double *cos);
void sincosl(long double x,long double *sin,long double *cos);
void sincosq(__float128 x,__float128 *sin,__float128 *cos);

inline int sign (val_type a)
{
 return ((a) >= 0 ? 1 : -1);
}

inline int sign_high (highprec a)
{
 return ((a) >= 0 ? 1 : -1);
}

/************************************************************************************/
/* 					   					    */
/* KeplerflowAll		         	  			            */
/*                                        					    */
/*										    */
/************************************************************************************/
void KeplerFlowAll (int neq,int keplerkop, solution *u,
                    val_type h, parameters *params)
{
#define BASE val_type
#define HIGH 0
#include <source_KeplerFlowAll.h>
#undef BASE
#undef HIGH
}


/************************************************************************************/
/* 					   					    */
/* KeplerflowAll2		         	  			            */
/* For stages initialization			             			    */
/*										    */
/************************************************************************************/
void KeplerFlowAll2 (int neq,int keplerkop, val_type *u,
                     val_type h, parameters *params)
{

#define BASE val_type
#define HIGH 0
#include <source_KeplerFlowAll2.h>
#undef BASE
#undef HIGH

}


/************************************************************************************/
/* 					   					    */
/* KeplerflowAll_high 		         	  			            */
/*                                        					    */
/*										    */
/************************************************************************************/
void KeplerFlowAll_high (int neq,int keplerkop, solution *u,
                         highprec h, parameters *params)
{
#define BASE highprec
#define HIGH 1
#include <source_KeplerFlowAll.h>
#undef BASE
#undef HIGH
}



/************************************************************************************/
/* 					   					    */
/* KeplerFlow			         	  			            */
/*                                        					    */
/*										    */
/************************************************************************************/

void KeplerFlow (val_type k, val_type *q,val_type *v,
                 val_type *dq,val_type *dv,val_type dt)
{
#define BASE val_type
#define KSQRT(x) SQRT(x)
#define HIGH 0
#include <source_KeplerFlow.h>
#undef BASE
#undef KSQRT
#undef HIGH
}


/************************************************************************************/
/* 					   					    */
/* KeplerFlow_high  		         	  			            */
/*                                        					    */
/*										    */
/************************************************************************************/

void KeplerFlow_high (highprec k, highprec *q,highprec *v,
                      highprec *dq,highprec *dv,highprec dt)
{
#define BASE highprec
#define KSQRT(x) SQRT_high(x)
#define HIGH 1
#include <source_KeplerFlow.h>
#undef BASE
#undef KSQRT
#undef HIGH
}



/************************************************************************************/
/* 					   					    */
/* KeplerSolveE			         	  			            */
/*                                        					    */
/*										    */
/************************************************************************************/
void KeplerSolveE (val_type r0, val_type eta, val_type zeta,  
                   val_type beta, val_type k, val_type t,
                   val_type *X, val_type *G, int *iter)
{
#define BASE val_type
#define HIGH 0
#define KSQRT(x) SQRT(x)
#define KSIN(x)  SIN(x)
#define KCOS(x)  COS(x)
#define KSINCOS(x,s,c)  SINCOS(x,s,c)
#define KFABS(x) FABS(x)
#define KPOW(x,y) POW(x,y)		
#define KTOL TOLKEPLER
#include <source_KeplerSolveE.h>
#undef BASE
#undef HIGH
#undef KSQRT 
#undef KSIN  
#undef KCOS 
#undef KSINCOS 
#undef KFABS 
#undef KPOW
#undef KTOL
}

/************************************************************************************/
/* 					   					    */
/* KeplerSolveE_high 		         	  			            */
/*                                        					    */
/*										    */
/************************************************************************************/
void KeplerSolveE_high (highprec r0, highprec eta, highprec zeta, 
                        highprec beta, highprec k, highprec t,
                        highprec *X, highprec *G, int *iter) 
{
#define LOW val_type
#define BASE highprec
#define HIGH 1
#define KSQRT(x) SQRT_high(x)
#define KSIN(x)  SIN_high(x)
#define KCOS(x)  COS_high(x)
#define KSINCOS(x,s,c)  SINCOS_high(x,s,c)
#define KFABS(x) FABS_high(x)
#define KPOW(x,y) POW_high(x,y)		
#define KTOL TOLKEPLER_high
#include <source_KeplerSolveE_high.h>
#undef LOW
#undef BASE
#undef HIGH
#undef KSQRT 
#undef KSIN  
#undef KCOS  
#undef KSINCOS
#undef KFABS
#undef KPOW
#undef KTOL

}



/************************************************************************************/
/* 					   					    */
/* KeplerFlowGFcn  					         	            */
/*                                        					    */
/*										    */
/************************************************************************************/

void  KeplerFlowGFcn     (void GFcn(), int neq, val_type t,
                         val_type *U, int keplerkop,
                         val_type *k, parameters *params, val_type dt,
                         val_type *G)
{
#define BASE val_type
#define HIGH 0
#include <source_KeplerFlowGFcn.h>
#undef BASE
#undef HIGH

}

void  KeplerFlowGFcn_high (void GFcn(), int neq, highprec t,
                           highprec *U, int keplerkop,
                           highprec *k, parameters_high *params, highprec dt,
                           highprec *G)
{
#define BASE highprec
#define HIGH 1
#include <source_KeplerFlowGFcn.h>
#undef BASE
#undef HIGH

}


/************************************************************************************/
/* 					   					    */
/* KeplerFlowGen 				 			            */
/*                                        					    */
/*										    */
/************************************************************************************/


void KeplerFlowGen    (val_type dt, val_type *Q, val_type *V, val_type k,
                       val_type *qq, val_type *vv, flowaux *aux)
{
#define BASE val_type
#define HIGH 0
#define KSQRT(x) SQRT(x)
#include <source_KeplerFlowGen.h>
#undef KSQRT
#undef BASE
#undef HIGH
}

void KeplerFlowGen_high (highprec dt, highprec *Q, highprec *V, highprec k,
                         highprec *qq, highprec *vv, flowaux_high *aux)
{
#define BASE highprec
#define HIGH 1
#define KSQRT(x) SQRT_high(x)
#include <source_KeplerFlowGen.h>
#undef KSQRT
#undef BASE
#undef HIGH
}


/************************************************************************************/
/* 					   					    */
/* KeplerFlowGFcnAux 							            */
/*                                        					    */
/*										    */
/************************************************************************************/

void KeplerFlowGFcnAux (val_type dt,val_type *Q, val_type *V, 
                        val_type *gq, val_type *gv,
                        flowaux *aux, val_type k, val_type *Gq, val_type *Gv)
{
#define BASE val_type
#define HIGH 0
#include <source_KeplerFlowGFcnAux.h>
#undef BASE
#undef HIGH
}

void KeplerFlowGFcnAux_high (highprec dt,highprec *Q, highprec *V, 
                        highprec *gq, highprec *gv,
                        flowaux_high *aux, highprec k, highprec *Gq, highprec *Gv)
{
#define BASE highprec
#define HIGH 1
#include <source_KeplerFlowGFcnAux.h>
#undef BASE
#undef HIGH
}


/************************************************************************************/
/* 					   					    */
/* ProjFun 			         	  			            */
/*                                        					    */
/*										    */
/************************************************************************************/


void ProjFun   (int neq, val_type t, val_type h, solution *w, parameters *params)
{
#define BASE val_type
#define HIGH 0
#include <source_ProjFun.h>
#undef BASE
#undef HIGH
}


/************************************************************************************/
/* 					   					    */
/* ProjFun_high 		         	  			            */
/*                                        					    */
/*										    */
/************************************************************************************/


void ProjFun_high (int neq, val_type t, val_type h, solution *w, parameters *params)
{
#define BASE highprec
#define HIGH 1
#include <source_ProjFun.h>
#undef BASE
#undef HIGH
}

/************************************************************************************/
/* 					   					    */
/* StartFun 			         	  			            */
/*                                        					    */
/*										    */
/************************************************************************************/


void StartFun  (int neq, val_type t, val_type h, solution *w, parameters *params)
{
#define BASE val_type
#define HIGH 0
#include <source_StartFun.h>
#undef BASE
#undef HIGH
}


void StartFun_high  (int neq, val_type t, val_type h, solution *w, parameters *params)
{
#define BASE highprec
#define HIGH 1
#include <source_StartFun.h>
#undef BASE
#undef HIGH
}

/******************************************************************************/
/*									      */
/*      OutputFun : 							      */
/*                   		                                              */
/*									      */
/******************************************************************************/
 
void OutputFun ( ode_sys *system, tcoeffs *method, val_type t, 
                  val_type h, solution *w, tcache_stat *cache_stat,
                  parameters *params, toptions *options, FILE *myfile)
{
#define BASE val_type
#define HIGH 0
#include <source_OutputFun.h>   
#undef BASE
#undef HIGH
}

void OutputFun_high ( ode_sys *system, tcoeffs *method, val_type t, 
                       val_type h, solution *w, tcache_stat *cache_stat,
                       parameters *params, toptions *options, FILE *myfile)
{
#define BASE highprec
#define HIGH 1
#include <source_OutputFun.h>   
#undef BASE
#undef HIGH
}



