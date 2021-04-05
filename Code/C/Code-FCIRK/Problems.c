/*------------------------------------------------------------------------------*/
/*										*/
/*    File: Problem.c								*/
/*										*/
/*        N-Body equation of the Solar System using Canonical Heliocentric      */
/*      Coordinates (point mass Newtonian model).                               */
/*                                                                              */
/*	    du/dt=k(u)+g(u), where						*/
/*             k(u): uncoupled Keplerians motions                               */
/*             g(u): interations among planets                                  */
/*										*/
/*	Functions:								*/
/*										*/
/*         HamNbody()= H_K+H_I                                                  */
/*                                                                              */
/*         H_K : Kepler Hamiltonian                                             */
/*             Ham_k: Hamiltonian                                               */
/*                                                                              */
/*         H_I : Interaction Hamiltonian                                        */
/*                                                                              */
/*             1-Problem :   Considering Earth+Moon one body                    */
/*                 RR1:      Interaction Hamiltonian                  		*/
/*                 DRR1:     Gradient of interaction Hamiltonian                */
/*                 KComp1:   K/Mu Computation                       		*/
/*         									*/
/*             2-Problem :   Considering Moon as separate body       		*/
/*                 RR2:      Interaction Hamiltonian                            */
/*                 DRR2:     Gradient of interaction Hamiltonian                */
/*                 KComp2:   K/Mu Computation                        	        */
/*                 Mysubstract: Implementation of substract                     */
/*                              to avoid cancellation                           */
/*         									*/
/*	   Transformed ODE system of N-Body problem		                */
/*             NbodyOde:     Transformed equations by                           */
/*                           Kepler flow's change variables                     */
/*             NbodyGFcn:    Gradient of interaction Hamiltonian                */
/*         									*/
/*------------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <def.h>
#include <Problems.h>
#include <quadmath.h>


/********************************************************************************/
/*										*/
/*          HamNbody                                                            */
/*             - H(Q,V)=H_K(Q,V)+H_I(Q,V)                                       */
/*             - Poincare's canonical Heliocentric coordinates                  */
/*                                                                              */
/********************************************************************************/


__float128 HamNbody (int neq,solution *u,parameters *params)
{

/* ---------- First initializations ------------------------------------*/

     int dim,nbody;
     dim=3;
     nbody=neq/(2*dim);

     Pkepler_sys *Pkepler;


/*------ declarations -------------------------------------------------*/


     __float128 H;


/* ----------- implementation  ----------------------------------------*/

     Pkepler=&params->Pkepler;
     H=Ham_K(nbody,u,params)+Pkepler->RR(neq,u,params);

#ifdef DEBUGODE3

     printf("Ham3\n");
     int i,n;
     int width = 46;
     char buf[128];
     __float128 aux;

     printf("Ham_K\n");
     aux=Ham_K(nbody,u,params);
     n = quadmath_snprintf(buf, sizeof buf, "%+-#*.20Qe", width, aux);
     if ((size_t) n < sizeof buf) printf("%s\n",buf);

     printf("PKepler\n");
     aux=Pkepler->RR(neq,u,params);
     n = quadmath_snprintf(buf, sizeof buf, "%+-#*.20Qe", width, aux);
     if ((size_t) n < sizeof buf) printf("%s\n",buf);


#endif 
     return(H);

}



/********************************************************************************/
/*										*/
/*              Ham_K                                                           */
/*                                                                              */
/********************************************************************************/


__float128 Ham_K (int nbody,solution *u, parameters *params)
{

/* ---------- First initializations ------------------------------------*/

     int dim;
     int neq;

/*------ declarations -------------------------------------------------*/

     __float128 *uu;
     int i,id,i1,i2;
     int nd;
     __float128 Vi2,Qi,H;


    Pkepler_sys *Pkepler;
    __float128 *k,*mu;

/* ----------- implementation  ----------------------------------------*/

     Pkepler=&params->Pkepler;
     k=Pkepler->Khigh;
     mu=Pkepler->Muhigh;

     dim=3;
     neq=nbody*2*dim;
     uu =(__float128 *)malloc(sizeof(__float128)*neq);

     for (i=0; i<neq; i++)   uu[i]=u->uu[i];

     nd=nbody*dim;
     H=0.0;

     for (i=0; i<nbody; i++)
     {
           i1=dim*i;
           i2=nd+i1;
           Vi2=0.0;
           Qi=0.0;

           for (id=0; id<dim; id++)
           {
                Vi2+=uu[i2+id]*uu[i2+id];
                Qi+= uu[i1+id]*uu[i1+id];
           }

           Qi=sqrtq(Qi);
           H+= (Vi2*mu[i]/2.)-(mu[i]*k[i]/Qi);
      }

     free(uu);

     return(H);

}



/***********************************************************************************/
/*										   */
/*	RR1: Interaction Hamiltonian	 					   */
/*										   */
/*										   */
/***********************************************************************************/


__float128 RR1 (int neq,solution *u,parameters *params)
{

/* ---------- First initializations ------------------------------------*/

     int dim,nbody;
     dim=3;
     nbody=neq/(2*dim);

/*------ declarations -------------------------------------------------*/

     __float128 uu[neq];
     int i,id,j,i1,i2,j1,j2;
     int nd;

     __float128 d,ViVj;
     __float128 T,U,R;

     __float128 *Gm;
     __float128 GmSun,Gma[nbody];

     __float128 *mu;
     Pkepler_sys *Pkepler;

/* ----------- implementation  ----------------------------------------*/

     Pkepler=&params->Pkepler;
     mu=Pkepler->Muhigh;

     nd=neq/2;

     for (i=0; i<neq; i++)	uu[i]=u->uu[i];
     Gm=params->rparhigh;

     GmSun=Gm[0];
     for (i=0; i<nbody; i++)    Gma[i]=Gm[i+1];

     T=0.0;
     U=0.0;

     for (i=0; i<nbody; i++)
     {
           i1=dim*i;
           i2=nd+i1;

           for (j=i+1; j<nbody; j++)
           {
               j1=j*dim;
               j2=nd+j1;
               d=0.0;
               for (id=0; id<dim; id++) d+=(uu[i1+id]-uu[j1+id])*(uu[i1+id]-uu[j1+id]);
               d=sqrtq(d);

               ViVj=0.0;
               for (id=0; id<dim; id++) ViVj+=uu[i2+id]*uu[j2+id];
               T+=(mu[i]*mu[j]*ViVj);
               U+=Gma[i]*Gma[j]/d;
           }
     }

     R=T/GmSun-U;
     return(R);

}


/*-----------------------------------------------------------------------------*/
/*									       */
/*                        DRR1			         		       */
/*									       */
/* ----------------------------------------------------------------------------*/

void DRR1 (int neq, val_type t,val_type *u,val_type *dR,parameters *params)
{
#define BASE val_type
#define HIGH 0
#define PSQRT(x) SQRT(x)
#include <source_DRR1.h>
#undef PSQRT
#undef BASE
#undef HIGH
}


void DRR1_high (int neq, highprec t,highprec *u,highprec *dR,parameters_high *params)
{
#define BASE highprec
#define HIGH 1
#define PSQRT(x) SQRT_high(x)
#include <source_DRR1.h>
#undef PSQRT
#undef BASE
#undef HIGH
}


/***********************************************************************************/
/*										   */
/*	KComp1: K/Mu parameters computation					   */
/*      N-Body problem 								   */
/*										   */
/***********************************************************************************/


void KComp1 (int nbody, ode_sys *system, Pkepler_sys *Pkepler)
{

/*------ Declarations --------------------------------------------------------*/

     int i;
     __float128 Gmsun, Gma[nbody];

     parameters *params;

/* ----- Implementation ------------------------------------------------------*/

     params=&system->params;

     Gmsun=params->rparhigh[0];

     
     for (i=0; i<nbody; i++)
     {
       Gma[i]=params->rparhigh[i+1];
       Pkepler->Muhigh[i]=(Gmsun*Gma[i])/(Gmsun+Gma[i]);
       Pkepler->K[i]=Pkepler->Khigh[i];
       Pkepler->Mu[i]=Pkepler->Muhigh[i];
       Pkepler->K2[i]=Pkepler->Khigh[i]-Pkepler->K[i];            
     }


#ifdef DEBUGKCOMP1

    printf("k values --------\n");
    double daux0,daux1,daux2;
    for (i=0; i<nbody; i++)
    {
        daux0=Pkepler->Khigh[i];
        daux1=Pkepler->K[i];
        daux2=Pkepler->K2[i];
        printf("nbody=%i, khigh=%.20lg, k=%.20lg, k2=%.20lg\n", 
                i, daux0,daux1, daux2);
    }  
#endif

    return;

}



/******************************************************************************/
/*									      */
/*	RR2: Interaction Hamiltonian					      */
/*	     Moon as separate body					      */
/*									      */
/******************************************************************************/


__float128 RR2 (int neq,solution *u,parameters *params)
{

/*------ Declarations --------------------------------------------------------*/

     int dim,nbody;
     dim=3;
     nbody=neq/(2*dim);

     __float128 uu[neq];
     int i,id,j,i1,i2,j1,j2;
     int iM,iM1,iE,iE1;
     int nd;

     __float128 da,db,dc,ViVj;
     __float128 T1,T2,T3,R;

     __float128 *Gm;
     __float128 GmSun,GmEM,Gma[nbody];

     __float128 qa[dim],qb[dim];

     Pkepler_sys *Pkepler;
     __float128 *mu;

/* ----- Implementation ------------------------------------------------------*/

     Pkepler=&params->Pkepler;
     mu=Pkepler->Muhigh;

     nd=neq/2;

     iM=nbody-1;  	     // Moon
     iM1=iM*dim;
     iE=nbody-2;             // Earth
     iE1=iE*dim;

     for (i=0; i<neq; i++) uu[i]=u->uu[i];


     Gm=params->rparhigh;
     GmSun=Gm[0];
     for (i=0; i<nbody; i++)  Gma[i]=Gm[i+1];
     GmEM=Gma[iM]/Gma[iE];

     T1=0;
     T2=0.;
     T3=0.;

     for (i=0; i<(nbody-1); i++)
     {
           i1=dim*i;
           i2=nd+i1;

           for (j=i+1; j<(nbody-1); j++)
           {
               j1=j*dim;
               j2=nd+j1;

               ViVj=0.0;
               for (id=0; id<dim; id++) ViVj+=uu[i2+id]*uu[j2+id];
               T1+=mu[i]*mu[j]*ViVj;
           }
     }


     for (i=0; i<(nbody-2); i++)
     {
           i1=dim*i;
           i2=nd+i1;

           for (j=i+1; j<(nbody-2); j++)
           {
               j1=j*dim;
               j2=nd+j1;
               da=0.;
               for (id=0; id<dim; id++) da+=(uu[i1+id]-uu[j1+id])*(uu[i1+id]-uu[j1+id]);
               da=sqrtq(da);
               T2+=Gma[i]*Gma[j]/da;
           }
     }


     /* Moon interaction */

     da=0.;
     db=0.;
     dc=0.;
     for (id=0; id<dim; id++)
     {
        qa[id]=uu[iE1+id]-GmEM*uu[iM1+id];
        qb[id]=uu[iE1+id]+uu[iM1+id];
        da+=qa[id]*qa[id];
        db+=qb[id]*qb[id];
        dc+=uu[iE1+id]*uu[iE1+id];
     }

     da=sqrtq(da);
     db=sqrtq(db);
     dc=sqrtq(dc);

     T3=GmSun*(Gma[iE]+Gma[iM])/dc-GmSun*Gma[iE]/da-GmSun*Gma[iM]/db;


     for (i=0; i<(nbody-2); i++)
     {
         i1=i*dim;
         i2=nd+i1;
         da=0.;
         db=0.;
         for (id=0; id<dim; id++)
         {
            qa[id]=(uu[i1+id]-uu[iE1+id]-uu[iM1+id]);
            qb[id]=(uu[i1+id]-uu[iE1+id]+GmEM*uu[iM1+id]);
            da+=qa[id]*qa[id];
            db+=qb[id]*qb[id];
         }

         da=sqrtq(da);
         db=sqrtq(db);

         T3-= Gma[i]*Gma[iM]/da+Gma[i]*Gma[iE]/db;

     }


     R=T1/GmSun-T2+T3;


#ifdef DEBUGRR2

     int n;
     int width = 46;
     char buf[128];

     printf("T1:");
     n = quadmath_snprintf(buf, sizeof buf, "%+-#*.4Qe", width, T1);
     if ((size_t) n < sizeof buf) printf("%s\n",buf);

     printf("T2:");
     n = quadmath_snprintf(buf, sizeof buf, "%+-#*.4Qe", width, T2);
     if ((size_t) n < sizeof buf) printf("%s\n",buf);

     printf("T3:");
     n = quadmath_snprintf(buf, sizeof buf, "%+-#*.4Qe", width, T3);
     if ((size_t) n < sizeof buf) printf("%s\n",buf);


#endif
     return(R);

}


/******************************************************************************/
/*									      */
/*	Mysubstract							      */
/*									      */
/******************************************************************************/


void Mysubstract (val_type *q, val_type *Dq, val_type *sub)
{

  /* Compute: (q+Dq)/||q+Dq||^3 - q/||q||^3                               */
  /*     as    (a^3-b^3)q+b^3*Dq,     a=1/||q+Dq||, b=1/||q||             */
  /*                                                                      */

#define BASE val_type
#define HIGH 0
#define PSQRT(x) SQRT(x)
#include <source_Mysubstract.h>
#undef PSQRT
#undef BASE
#undef HIGH

}



void Mysubstract_high (highprec *q, highprec *Dq, highprec *sub)
{

  /* Compute: (q+Dq)/||q+Dq||^3 - q/||q||^3                               */
  /*     as    (a^3-b^3)q+b^3*Dq,     a=1/||q+Dq||, b=1/||q||             */
  /*                                                                      */

#define BASE highprec
#define HIGH 1
#define PSQRT(x) SQRT_high(x)
#include <source_Mysubstract.h>
#undef PSQRT
#undef BASE
#undef HIGH

}


/******************************************************************************/
/*									      */
/*	DRR2:  Gradient of interaction Hamiltonian			      */
/*	       Moon as separate body					      */
/*									      */
/******************************************************************************/
void DRR2 (int neq, val_type t,val_type *u,val_type *dR,parameters *params)
{

#define BASE val_type
#define HIGH 0
#define PSQRT(x) SQRT(x)
#include <source_DRR2.h>
#undef PSQRT
#undef BASE
#undef HIGH

}

void DRR2_high (int neq, highprec t,highprec *u,highprec *dR,parameters_high *params)
{

#define BASE highprec
#define HIGH 1
#define PSQRT(x) SQRT_high(x)
#include <source_DRR2.h>
#undef PSQRT
#undef BASE
#undef HIGH

}


/***********************************************************************************/
/*										   */
/*	KComp2: K/Mu parameters computation					   */
/*      									   */
/*										   */
/***********************************************************************************/


void KComp2 (int nbody, ode_sys *system, Pkepler_sys *Pkepler)
{

/*------ Declarations --------------------------------------------------------*/

     int i,iE,iM;
     __float128 Gmsun, Gma[nbody];

     parameters *params;

/* ----- Implementation ------------------------------------------------------*/

     params=&system->params;

     iE=nbody-2;
     iM=nbody-1;

     Gmsun=params->rparhigh[0];
     
     for (i=0; i<nbody; i++)
     {
       Gma[i]=params->rparhigh[i+1];
       Pkepler->Muhigh[i]=(Gmsun*Gma[i])/(Gmsun+Gma[i]);
       Pkepler->K[i]=Pkepler->Khigh[i];
       Pkepler->Mu[i]=Pkepler->Muhigh[i];
       Pkepler->K2[i]=Pkepler->Khigh[i]-Pkepler->K[i];            
     }


     Pkepler->Muhigh[iE]=Gmsun*(Gma[iE]+Gma[iM])/(Gmsun+Gma[iE]+Gma[iM]);
     Pkepler->Muhigh[iM]=Gma[iM]*(Gma[iE]+Gma[iM])/Gma[iE];
     Pkepler->K[iE]=Pkepler->Khigh[iE];
     Pkepler->Mu[iE]=Pkepler->Muhigh[iE];
     Pkepler->K[iM]=Pkepler->Khigh[iM];
     Pkepler->Mu[iM]=Pkepler->Muhigh[iM];
     Pkepler->K2[iE]=Pkepler->Khigh[iE]-Pkepler->K[iE];            
     Pkepler->K2[iM]=Pkepler->Khigh[iM]-Pkepler->K[iM];            

#ifdef DEBUGKCOMP2
    printf("k values --------\n");
    double daux0,daux1,daux2;
    for (i=0; i<nbody; i++)
    {
        daux0=Pkepler->Khigh[i];
        daux1=Pkepler->K[i];
        daux2=Pkepler->K2[i];
        printf("nbody=%i, khigh=%.20lg, k=%.20lg, k2=%.20lg\n", i, daux0,daux1, daux2);
    }  
#endif


    return;

}




/***********************************************************************************/
/*										   */
/*	NbodyGFcn 								   */
/*										   */
/***********************************************************************************/


void NbodyGFcn (int neq, val_type t,val_type *u,val_type *dR,parameters *params)
{

#define BASE val_type
#define HIGH 0
#include <source_NbodyGFcn.h>
#undef BASE
#undef HIGH

}


void NbodyGFcn_high (int neq, highprec t,highprec *u,highprec *dR,parameters_high *params)
{

#define BASE highprec
#define HIGH 1
#include <source_NbodyGFcn.h>
#undef BASE
#undef HIGH

}

/***********************************************************************************/
/*										   */
/*	NbodyOde: transformed equations by Kepler flow's change variables	   */
/*										   */
/***********************************************************************************/

void NbodyOde (int neq, val_type t,val_type ttau,
               val_type *u,val_type *f,parameters *params)
{

#define BASE val_type
#define HIGH 0
#define PSQRT(x) SQRT(x)
#include <source_NbodyOde.h>
#undef PSQRT
#undef BASE
#undef HIGH

}

void NbodyOde_high (int neq, highprec t,highprec ttau,
                    highprec *u,highprec *f,parameters_high *params)
{

#define BASE highprec
#define HIGH 1
#define PSQRT(x) SQRT_high(x)
#include <source_NbodyOde.h>
#undef PSQRT
#undef BASE
#undef HIGH

}


