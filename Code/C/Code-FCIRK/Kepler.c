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

#include <source_KeplerFlowAll2.h>

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

     int dim=3;
     int nd=neq/2;

     int j,j1,j2;

     val_type u[neq],g[neq];
     flowaux aux[keplerkop];


     for (j=0; j<neq; j++) u[j]=U[j];

     for (j=0; j<keplerkop; j++)
     {
          j1=j*dim;
          j2=j1+nd;
          KeplerFlowGen (dt, &U[j1], &U[j2], k[j], &u[j1], &u[j2], &aux[j]);

     }

     GFcn (neq, t, u, g, params);
     for (j=0; j<neq; j++) G[j]=g[j];

     for (j=0; j<keplerkop; j++)
     {
          j1=j*dim;
          j2=j1+nd;
          KeplerFlowGFcnAux (dt,&U[j1],&U[j2],&g[j1],&g[j2],&aux[j],k[j],&G[j1],&G[j2]);

     }


#    ifdef DEBUGKEPLERFLOWGFCN

   printf("KeplerFlowGFcn, neq=%i\n",neq);
   
   int i;
   val_type norma;

   printf("dt=%lg\n",dt);

   norma=0.;
   for (i=0; i<neq-1; i++) norma+=G[i]*G[i];
   printf("Norm(f)=%lg\n",sqrt(norma));
   printf("\n");

#    endif

   return;

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

     int id;
     int dim=3;   


#ifdef XDEBUG         
     printf("KeplerFlowGen.....");    
               
     double q0,q1,q2;
     double v0,v1,v2;
     q0=Q[0]; q1=Q[1]; q2=Q[2];
     v0=V[0]; v1=V[1]; v2=V[2];
     printf("KeplerFlowGen   Q=(%lg,%lg,%lg), V=(%lg,%lg,%lg) !!!\n",q0,q1,q2,v0,v1,v2);

#endif


     if (k==0)
     {

        for (id=0; id<dim; id++)
        {
          qq[id]=Q[id]+dt*V[id];
          vv[id]=V[id];
        }

     }
     else
     {    
          
        aux->r0=SQRT(Q[0]*Q[0]+Q[1]*Q[1]+Q[2]*Q[2]);
        aux->eta=Q[0]*V[0]+Q[1]*V[1]+Q[2]*V[2];
        aux->alpha=k/(aux->r0);
        aux->beta=2*(aux->alpha)-(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
        aux->gamma=k/(aux->beta);
        aux->zeta=k-(aux->beta)*(aux->r0);

        KeplerSolveE(aux->r0,aux->gamma,aux->eta,aux->beta,k,dt,&aux->X,aux->GG);	
        aux->r= (aux->r0)+(aux->eta)*aux->GG[1]+(aux->zeta)*aux->GG[2];	        
        aux->rinv=1/(aux->r);
     
        aux->b[0]=-(aux->alpha)*aux->GG[2];
        aux->b[1]=(aux->r0)*aux->GG[1]+(aux->eta)*aux->GG[2];
        aux->b[2]=-(aux->alpha)*aux->GG[1]*(aux->rinv);
        aux->b[3]=-k*aux->GG[2]*(aux->rinv);

        for (id=0; id<dim; id++)
        {
           qq[id]=(1+aux->b[0])*Q[id]+aux->b[1]*V[id];
           vv[id]=(1+aux->b[3])*V[id]+aux->b[2]*Q[id];
        }

     }


#    ifdef MDEBUG

     int i;
     double xx;

     printf("KeplerFlowGen*******************\n");


     printf("dt=%lg\n",dt);

     printf("r0=%lg,eta=%lg,alpha=%lg,beta=%lg,gamma=%lg,zeta=%lg\n",
             aux->r0,aux->eta,aux->alpha,aux->beta,aux->gamma,aux->zeta);
     printf("X=%lg,G0=%lg, G1=%lg, G2=%lg, r=%lg, rinv=%lg\n",
             aux->X,aux->GG[0],aux->GG[1],aux->GG[2],aux->r,aux->rinv);
     printf("b11=%lg, b12=%lg, b21=%lg, b22=%lg\n", 
             aux->b[0],aux->b[1],aux->b[2],aux->b[3]);

     printf("q,v=");
     for (i=0; i<dim; i++)
     {  xx=qq[i]; 
        printf("%lg,",xx);
     }
     printf("\n");
     for (i=0; i<dim; i++)
     {
        xx=vv[i];
        printf("%lg,",xx);
     }
     printf("\n");

#    endif

     return;

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

     int dim=3;

     int id;

     val_type r0inv,betainv;
     val_type lambda[16];
     val_type r0kg,G1b,G2b;
     val_type Rqv[2*dim],RQV[2*dim];


     if (k==0)			
     {

        for (id=0; id<dim; id++)
        {
          Gq[id]=gq[id]-dt*gv[id];		
          Gv[id]=gv[id];	
		
        }


     }
     else
     {


        for (id=0; id<dim; id++)
        {
           Rqv[id]=-gv[id];        // (Rq1,Rq2,Rq3)<---(-gv1,-gv2.-gv3)
           Rqv[dim+id]=gq[id];     // (Rv1,Rv2,Rv3)<---( gq1, gq2, gq3)

        }

    
        betainv=1/aux->beta;
        r0inv=1/aux->r0;
        G1b=0.5*(aux->X*aux->GG[0]-aux->GG[1])*betainv;    
        G2b=(0.5*aux->X*aux->GG[1]-aux->GG[2])*betainv;    

        lambda[15]=(-V[0]*Rqv[dim]-V[1]*Rqv[dim+1]-V[2]*Rqv[dim+2])*aux->rinv;
        lambda[14]=(-Q[0]*Rqv[dim]-Q[1]*Rqv[dim+1]-Q[2]*Rqv[dim+2])*aux->rinv;
        lambda[13]=(-V[0]*Rqv[0]-V[1]*Rqv[1]-V[2]*Rqv[2]);
        lambda[12]=-Q[0]*Rqv[0]-Q[1]*Rqv[1]-Q[2]*Rqv[2];
        lambda[11]=aux->b[2]*lambda[14]+aux->b[3]*lambda[15];
        lambda[10]=-aux->zeta*lambda[11]-aux->alpha*lambda[12]+aux->eta*lambda[13]-k*lambda[15];
        lambda[9]=-aux->eta*lambda[11]+aux->r0*lambda[13]-aux->alpha*lambda[14];
        r0kg=aux->r0-aux->gamma;
        lambda[8]=(lambda[9]*aux->GG[0]+lambda[10]*aux->GG[1])/(aux->gamma+r0kg*aux->GG[0]+aux->eta*aux->GG[1]);
        lambda[7]=-aux->GG[1]*lambda[8];
        lambda[6]=-aux->GG[2]*lambda[11];
        lambda[5]=(-lambda[7]-aux->X*lambda[8])*betainv;
        lambda[4]=2*(-aux->gamma*lambda[5]-aux->r0*lambda[6]+
                  (lambda[9]-r0kg*lambda[8])*G1b+(lambda[10]-aux->eta*lambda[8])*G2b);
        lambda[3]=(lambda[4]-aux->GG[2]*lambda[12]-aux->GG[1]*lambda[14])*r0inv;
        lambda[2]=(-aux->GG[1]*lambda[11]+aux->GG[2]*(lambda[13]-lambda[8]));
        lambda[1]=(-aux->alpha*lambda[3]-aux->beta*lambda[6]+lambda[7]-lambda[11]+aux->GG[1]*lambda[13])*r0inv;
 
   
        RQV[0]=-Q[0]*lambda[1]-V[0]*lambda[2]+(1+aux->b[0])*Rqv[0]+aux->b[2]*Rqv[dim];
        RQV[1]=-Q[1]*lambda[1]-V[1]*lambda[2]+(1+aux->b[0])*Rqv[1]+aux->b[2]*Rqv[dim+1];
        RQV[2]=-Q[2]*lambda[1]-V[2]*lambda[2]+(1+aux->b[0])*Rqv[2]+aux->b[2]*Rqv[dim+2];
        RQV[dim]=-Q[0]*lambda[2]+V[0]*lambda[4]+aux->b[1]*Rqv[0]+(1+aux->b[3])*Rqv[dim];
        RQV[dim+1]=-Q[1]*lambda[2]+V[1]*lambda[4]+aux->b[1]*Rqv[1]+(1+aux->b[3])*Rqv[dim+1];
        RQV[dim+2]=-Q[2]*lambda[2]+V[2]*lambda[4]+aux->b[1]*Rqv[2]+(1+aux->b[3])*Rqv[dim+2];

        for (id=0; id<dim; id++)
        {
             Gq[id]=RQV[dim+id];	 // (dQ1,dQ2,dQ3)<---(RV1,RV2.RV3)
             Gv[id]=-RQV[id];	 	 // (dV1,dV2,dV3)<---(-RQ1,-RQ2,-RQ3)
        }

#ifdef MDEBUG
        if (isnan(Gq[0]))
        {
           printf ("NAN KeplerFlowGFcnAux......\n");
           double dtx,Q1x, Q2x, Q3x;
           double V1x, V2x, V3x; 
           double gq1x,gq2x,gq3x;
           double gv1x,gv2x,gv3x;
           double kx;

           dtx=dt;
           Q1x=Q[0]; Q2x=Q[1]; Q3x=Q[2];
           V1x=V[0]; V2x=V[1]; V3x=V[2];
           gq1x=gq[0]; gq2x=gq[1]; gq3x=gq[2];
           gv1x=gv[0]; gv2x=gv[1]; gv3x=gv[2];
           kx=k;

           printf("dt=%lg, Q=%lg%lg,%lg, V=%lg,%lg,%lg\n", dtx, Q1x,Q2x,Q3x,V1x,V2x,V3x);
           printf("gq=%lg,%lg,%lg, gv=%lg,%lg,%lg, k=%lg\n",gq1x,gq2x,gq3x,gv1x,gv2x,gv3x,kx);  
        }
#endif
     }


     return;

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
 
void OutputFun ( ode_sys *system, gauss_method *method, val_type t, 
                  val_type h, solution *w, solver_stat *thestatptr,
                  parameters *params, toptions *options, FILE *myfile)
{
#define BASE6 val_type
#define HIGH6 0
#include <source_OutputFun.h>   
#undef BASE6
#undef HIGH6
}

void OutputFun_high ( ode_sys *system, gauss_method *method, val_type t, 
                       val_type h, solution *w, solver_stat *thestatptr,
                       parameters *params, toptions *options, FILE *myfile)
{
#define BASE6 highprec
#define HIGH6 1
#include <source_OutputFun.h>   
#undef BASE6
#undef HIGH6
}



