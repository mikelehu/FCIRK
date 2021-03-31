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
/*         KeplerFlowMIXED_high							*/
/*	   KeplerSolveE								*/
/*	   KeplerSolveE_high							*/
/*	   KeplerSolveITER_high							*/
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
/*	   source_KeplerFlowMIXED.c					        */
/*	   source_KeplerSolveE.c						*/
/*	   source_KeplerSolveITER.c						*/
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
#define BASE2 val_type
#define HIGH2 0
#include <source_KeplerFlowAll.h>
#undef BASE2
#undef HIGH2
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

#define BASE2B val_type
#define HIGH2B 0
#include <source_KeplerFlowAll2.h>
#undef BASE2B
#undef HIGH2B

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
#define BASE2 highprec
#define HIGH2 1
#include <source_KeplerFlowAll.h>
#undef BASE2
#undef HIGH2
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
#define BASE3 val_type
#define KSQRT(x) SQRT(x)
#define HIGH3 0
#include <source_KeplerFlow.h>
#undef BASE3
#undef KSQRT
#undef HIGH3
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
#define BASE3 highprec
#define KSQRT(x) SQRT_high(x)
#define HIGH3 1
#include <source_KeplerFlow.h>
#undef BASE3
#undef KSQRT
#undef HIGH3
}


/************************************************************************************/
/* 					   					    */
/* KeplerFlowMIXED_high  		       	  			            */
/*                                        					    */
/*										    */
/************************************************************************************/

void KeplerFlowMIXED_high (highprec k, highprec *q,highprec *v,
                         highprec *dq,highprec *dv,highprec dt)
{
#define BASE3 highprec
#define KSQRT(x) SQRT_high(x)
#define HIGH3 1
#include <source_KeplerFlowMIXED.h>
#undef BASE3
#undef KSQRT
#undef HIGH3
}

/************************************************************************************/
/* 					   					    */
/* KeplerSolveE			         	  			            */
/*                                        					    */
/*										    */
/************************************************************************************/
void KeplerSolveE (val_type r0, val_type gamma0, val_type eta0, 
                   val_type beta, val_type k, val_type t,
                   val_type *X, val_type *G)
{
#define BASE4 val_type
#define HIGH4 0
#define KSQRT(x) SQRT(x)
#define KSIN(x)  SIN(x)
#define KCOS(x)  COS(x)
#define KSINCOS(x,s,c)  SINCOS(x,s,c)
#define KFABS(x) FABS(x)
#define KPOW(x,y) POW(x,y)		
#define KTOL TOLKEPLER
#include <source_KeplerSolveE.h>
#undef BASE4
#undef HIGH4
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
void KeplerSolveE_high (highprec r0, highprec gamma0, highprec eta0, 
                        highprec beta, highprec k, highprec t,
                        highprec *X, highprec *G) 
{
#define BASE4 highprec
#define HIGH4 1
#define KSQRT(x) SQRT_high(x)
#define KSIN(x)  SIN_high(x)
#define KCOS(x)  COS_high(x)
#define KSINCOS(x,s,c)  SINCOS_high(x,s,c)
#define KFABS(x) FABS_high(x)
#define KPOW(x,y) POW_high(x,y)		
#define KTOL TOLKEPLER_high
#include <source_KeplerSolveE.h>
#undef BASE4
#undef HIGH4
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
/* KeplerSolveITER_high 		      	  			            */
/*                                        					    */
/*										    */
/************************************************************************************/


void KeplerSolveITER_high (highprec r0, highprec gamma, highprec eta, 
                            highprec beta, highprec k, highprec dt,
                            highprec *X, highprec *G) 
{
#define BASE4 highprec
#define HIGH4 1
#define KSQRT(x) SQRT_high(x)
#define KSIN(x)  SIN_high(x)
#define KCOS(x)  COS_high(x)
#define KFABS(x) FABS_high(x)
#define KTOL TOLKEPLER_high
#include <source_KeplerSolveITER.h>
#undef BASE4
#undef HIGH4
#undef KSQRT 
#undef KSIN  
#undef KCOS  
#undef KFABS
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
#define BASEFWGF val_type
#define HIGHFWGF 0
#include <source_KeplerFlowGFcn.h>
#undef BASEFWGF
#undef HIGHFWGF

}

void  KeplerFlowGFcn_high (void GFcn(), int neq, highprec t,
                           highprec *U, int keplerkop,
                           highprec *k, parameters_high *params, highprec dt,
                           highprec *G)
{
#define BASEFWGF highprec
#define HIGHFWGF 1
#include <source_KeplerFlowGFcn.h>
#undef BASEFWGF
#undef HIGHFWGF

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
#define BASEFGEN val_type
#define HIGHFGEN 0
#define KSQRT(x) SQRT(x)
#include <source_KeplerFlowGen.h>
#undef KSQRT
#undef BASEFGEN
#undef HIGHFGEN
}

void KeplerFlowGen_high (highprec dt, highprec *Q, highprec *V, highprec k,
                         highprec *qq, highprec *vv, flowaux_high *aux)
{
#define BASEFGEN highprec
#define HIGHFGEN 1
#define KSQRT(x) SQRT_high(x)
#include <source_KeplerFlowGen.h>
#undef KSQRT
#undef BASEFGEN
#undef HIGHFGEN
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
#define BASEFGAux val_type
#define HIGHFGAux 0
#include <source_KeplerFlowGFcnAux.h>
#undef BASEFGAux
#undef HIGHFGAux
}

void KeplerFlowGFcnAux_high (highprec dt,highprec *Q, highprec *V, 
                        highprec *gq, highprec *gv,
                        flowaux_high *aux, highprec k, highprec *Gq, highprec *Gv)
{
#define BASEFGAux highprec
#define HIGHFGAux 1
#include <source_KeplerFlowGFcnAux.h>
#undef BASEFGAux
#undef HIGHFGAux
}


/************************************************************************************/
/* 					   					    */
/* ProjFun 			         	  			            */
/*                                        					    */
/*										    */
/************************************************************************************/


void ProjFun   (int neq, val_type t, val_type h, solution *w, parameters *params)
{
#define BASE1 val_type
#define HIGH1 0
#include <source_ProjFun.h>
#undef BASE1
#undef HIGH1
}


/************************************************************************************/
/* 					   					    */
/* ProjFun_high 		         	  			            */
/*                                        					    */
/*										    */
/************************************************************************************/


void ProjFun_high (int neq, val_type t, val_type h, solution *w, parameters *params)
{
#define BASE1 highprec
#define HIGH1 1
#include <source_ProjFun.h>
#undef BASE1
#undef HIGH1
}

/************************************************************************************/
/* 					   					    */
/* StartFun 			         	  			            */
/*                                        					    */
/*										    */
/************************************************************************************/


void StartFun  (int neq, val_type t, val_type h, solution *w, parameters *params)
{
#define BASE5 val_type
#define HIGH5 0
#include <source_StartFun.h>
#undef BASE5
#undef HIGH5
}


void StartFun_high  (int neq, val_type t, val_type h, solution *w, parameters *params)
{
#define BASE5 highprec
#define HIGH5 1
#include <source_StartFun.h>
#undef BASE5
#undef HIGH5
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
#define BASE6 val_type
#define HIGH6 0
#include <source_OutputFun.h>   
#undef BASE6
#undef HIGH6
}

void OutputFun_high ( ode_sys *system, tcoeffs *method, val_type t, 
                       val_type h, solution *w, tcache_stat *cache_stat,
                       parameters *params, toptions *options, FILE *myfile)
{
#define BASE6 highprec
#define HIGH6 1
#include <source_OutputFun.h>   
#undef BASE6
#undef HIGH6
}



