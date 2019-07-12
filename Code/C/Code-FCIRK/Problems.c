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
/*         Hamkepler, Odekepler: Hamiltonian and Odefun of Keplerian part	*/
/*										*/
/*         Ode3()= Generic odefun : Odekepler+DRR  (IRKFP, CO1035)              */
/*         Ham3()= Generic Hamiltonian : Hkepler+RR                             */
/*										*/
/*	   NbodyOde: integration based on Kepler flow's change variables	*/
/*         NbodyGFcn                                                            */
/*                                                                              */
/*         1-Problem : Earth+Moon one body                                      */
/*                                                                              */
/*             RR1,DRR1: 10-Body perturbations (Ham and Gradiant)     		*/
/*             KComp1:   (K/Mu Computation)                            		*/
/*         									*/
/*         2-Problem : Moon as separate body					*/
/*         									*/
/*              RR2,DRR2: 11-Body perturbations (Ham and Gradiant)   		*/
/*              mysubstract                                                     */
/*              KComp2:   (K/Mu Computation)                           	        */
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
/*              HamKepler                                                       */
/*                                                                              */
/********************************************************************************/


__float128 HamKepler (int nbody,solution *u, parameters *params)
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


/********************************************************************************/
/*										*/
/*              OdeKepler                                                       */
/*                                                                              */
/********************************************************************************/


void OdeKepler (int neq, val_type t,val_type *u,val_type *f,parameters *params)
{

/* ---------- First initializations ------------------------------------*/

     int dim,nbody,nd;
     dim=3;
     nbody=neq/(2*dim);
     nd=neq/2;

/*------ declarations -------------------------------------------------*/

     int i,id,i1,i2;
     val_type Q3;
     val_type *k;

     Pkepler_sys *Pkepler;

/* ----------- implementation  ---------------------------------------*/


     Pkepler=&params->Pkepler;

     k=Pkepler->K;

     for (i=0; i<nbody; i++)
     {

           i1=i*dim;
           i2=nd+i1;

           Q3=0.;
           for (id=0; id<dim; id++) Q3+=u[i1+id]*u[i1+id];
           Q3=SQRT(Q3)*Q3;

           for (id=0; id<dim; id++)
           {
                f[i1+id]=u[i2+id];
                f[i2+id]=-(k[i]*u[i1+id])/Q3;
           }

      }


#ifdef DEBUGODEKEPLER
    printf("\nOdekepler\n");
    for (i=0; i<neq;i++) printf("%lg,", f[i]);
    printf("\n");
#endif

    return ;

}


/*------------------------------------------------------------------------------*/
/*										*/
/*        Problem-3:	 							*/
/*	    Ode3()= Generic odefun : Odekepler+DRR  (IRKFP, CO1035)		*/
/*	    Ham3()=: Generic Hamiltonian : Hankepler+RR				*/
/*										*/
/*------------------------------------------------------------------------------*/


void Ode3 (int neq, val_type t,val_type *u,val_type *f,parameters *params)

{
/* ---------- First initializations ------------------------------------*/

     int dim,nbody,nd;
     dim=3;
     nbody=neq/(2*dim);
     nd=neq/2;

/*------ declarations -------------------------------------------------*/

     int i,id,i1,i2;
     val_type dR[neq];

     Pkepler_sys *Pkepler;


/* ----------- implementation  ---------------------------------------*/

     Pkepler=&params->Pkepler;

     Pkepler->DRR (neq, t,u,dR,params);
     Pkepler->DK (neq,t,u,f,params);

     for (i=0; i<nbody; i++)
     {
          i1=i*dim;
          i2=nd+i1;
          for (id=0; id<dim; id++)
          {
               f[i1+id]+=dR[i1+id];
               f[i2+id]+=dR[i2+id];
          }
     }


#ifdef DEBUGODE3

     printf("Ode3\n");
     int n,k;
     int width = 20;  
     char buf[128];
     __float128 sum;
 

      printf("fq***\n");
      k=1;
      for (i=neq/2; i<neq; i++)
      {
           printf("%.20lg,",f[i]);
           if (k==3)
           {
              printf("\n");
              k=1;}
           else
           {  k++;}
       }

#endif

    return ;

}


__float128 Ham3 (int neq,solution *u,parameters *params)
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

     H=HamKepler(nbody,u,params)+Pkepler->RR(neq,u,params);


#ifdef DEBUGODE3

     printf("Ham3\n");
     int i,n;
     int width = 46;
     char buf[128];
     __float128 aux;

     printf("HamKepler\n");
     aux=HamKepler(nbody,u,params);
     n = quadmath_snprintf(buf, sizeof buf, "%+-#*.20Qe", width, aux);
     if ((size_t) n < sizeof buf) printf("%s\n",buf);

     printf("PKepler\n");
     aux=Pkepler->RR(neq,u,params);
     n = quadmath_snprintf(buf, sizeof buf, "%+-#*.20Qe", width, aux);
     if ((size_t) n < sizeof buf) printf("%s\n",buf);


#endif 
     return(H);

}



/***********************************************************************************/
/*										   */
/*	NbodyGFcn 								   */
/*										   */
/***********************************************************************************/


void NbodyGFcn (int neq, val_type t,val_type *u,val_type *dR,parameters *params)
{


/*------ Declarations --------------------------------------------------------*/

     Pkepler_sys *Pkepler;


/* ----- Implementation ------------------------------------------------------*/


     Pkepler=&params->Pkepler;
     Pkepler->DRR(neq,t,u,dR,params);
 

     return;

}


/***********************************************************************************/
/*										   */
/*	NbodyOde: transformed equations by Kepler flow's change variables	   */
/*										   */
/***********************************************************************************/

void NbodyOde (int neq, val_type t,val_type ttau,val_type *u,val_type *f,parameters *params)
{


/*------ Declarations --------------------------------------------------------*/

   val_type *k;
   Pkepler_sys *Pkepler;

/* ----- Implementation ------------------------------------------------------*/

   Pkepler=&params->Pkepler;
   k=Pkepler->K;   

   KeplerFlowGFcn (NbodyGFcn,neq,t,u,Pkepler->keplerkop,k,params,ttau,f);   


#    ifdef DEBUGNBODYODE

   printf("\nNbodyODe, neq=%i,keplerkop=%i\n",neq,Pkepler->keplerkop);

   int i;
   double aux;
   aux=ttau;
   printf("ttau=%lg\n",aux);

   aux=0.;
   for (i=0; i<neq-1; i++) aux+=ww[i]*ww[i];
   printf("Norm(w)=%lg\n",sqrt(aux));
   aux=0;
   for (i=0; i<neq-1; i++) aux+=f[i]*f[i];
   printf("Norm(f)=%lg\n",sqrt(aux));
   printf("\n");

#  endif

   return;

};



/***********************************************************************************/
/*										   */
/*	RR1: Hamiltonian of the perturbation 					   */
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
//     uu=u->uu;

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

/*------ Declarations --------------------------------------------------------*/

     int dim,nbody;
     dim=3;
     nbody=neq/(2*dim);

     int i,id,i1,i2,j,j1,j2;
     int nd;
     val_type d3,da,qij[dim];
     val_type *Gm,GmSun,Gma[nbody];

     Pkepler_sys *Pkepler;
     val_type *mu,*k2;

/* ----- Implementation ------------------------------------------------------*/

     nd=neq/2;
     Gm=params->rpar;
     Pkepler=&params->Pkepler;
     mu=Pkepler->Mu;
     k2=Pkepler->K2;

     GmSun=Gm[0];
     for (i=0; i<nbody; i++) {Gma[i]=Gm[i+1];}

     for (i=0; i<nbody; i++)
     {
           i1=i*dim;
           i2=nd+i1;
           for (id=0; id<dim; id++)
           {
                 dR[i1+id]=0.;
                 dR[i2+id]=0.;
           }
     }

     for (i=0; i<nbody; i++)
     {
          i1=i*dim;
          i2=nd+i1;

          for (j=i+1; j<nbody; j++)
          {
                j1=j*dim;
                j2=nd+j1;
                d3=0.;
                for (id=0; id<dim; id++)
                {
                     qij[id]=(u[i1+id]-u[j1+id]);
                     d3+=qij[id]*qij[id];
                }

                d3=SQRT(d3)*d3;

                for (id=0; id<dim; id++)
                {
                     dR[i1+id]+=mu[j]*u[j2+id];
                     dR[j1+id]+=mu[i]*u[i2+id];
                     dR[i2+id]-=Gma[j]*qij[id]/d3;
                     dR[j2+id]+=Gma[i]*qij[id]/d3;
                }
           }

          for (id=0; id<dim; id++)
          {
                dR[i1+id]*=1/GmSun;
                dR[i2+id]=(Gma[i]*dR[i2+id])/mu[i];
          }

     }


     for (i=0; i<nbody; i++)
     {

         if (k2[i]!=0)
         {
              i1=i*dim;
              i2=nd+i1;

              da=0.;
              for (id=0; id<dim; id++) da+=u[i1+id]*u[i1+id];
              da=SQRT(da)*da;
                 
              for (id=0; id<dim; id++)
                   dR[i2+id]-=k2[i]*u[i1+id]/da;
 
         }

     }
        


     return ;

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

    printf("k-ren balioak --------\n");
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
/*	RR2: Hamiltonian of the perturbation (N-Body)			      */
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


     /* Perturbation of the moon */

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
/*	mysubstract							      */
/*									      */
/******************************************************************************/


void mysubstract (val_type *q, val_type *Dq, val_type *sub)
{

  /* Compute: (q+Dq)/||q+Dq||^3 - q/||q||^3                               */
  /*     as    (a^3-b^3)q+b^3*Dq,     a=1/||q+Dq||, b=1/||q||             */
  /*                                                                      */

/*------ Declarations --------------------------------------------------------*/

    int dim=3;

    int id;
    val_type a3;

    val_type a,b;
    val_type l1,l2,l3,l4,l5,l6,l7,l8,l9;

/* ----- Implementation ------------------------------------------------------*/

    l1=0.;
    l2=0.;
    l7=0.;

    for (id=0; id<dim; id++)
    {
      l1+=q[id]*q[id];
      l2+=(q[id]+Dq[id])*(q[id]+Dq[id]);
      l7-=(2*q[id]+Dq[id])*Dq[id];
    }

    l3=SQRT(l1);
    l4=SQRT(l2);
    b=1./l3;
    a=1./l4;
    a3=1./(l4*l2);

    l5=a+b;
    l6=l5*a+1./l1;   //a^2+ab+b^2

    l8=1./(l3+l4);
    l9=l7*l8*a*b*l6;

    for (id=0; id<dim; id++) sub[id]=l9*q[id]+a3*Dq[id];


#ifdef DEBUGRR2

    for (id=0; id<dim; id++) printf("%.20lg,",sub[id]);
    printf("\n");

#endif

    return;

}



/******************************************************************************/
/*									      */
/*	DRR2:  Hamiltonian of the perturbation (N-Body)			      */
/*	       With the Moon						      */
/*									      */
/******************************************************************************/
void DRR2 (int neq, val_type t,val_type *u,val_type *dR,parameters *params)
{

/* ---------- First initializations ------------------------------------*/

     int dim,nbody;
     dim=3;
     nbody=neq/(2*dim);

/*------ declarations -------------------------------------------------*/

     int i,id,i1,i2,j,j1,j2;
     int iM,iM1,iM2,iE,iE1,iE2;
     int nd;
     val_type aux1,aux2;
     val_type qa[dim],qb[dim];
     val_type qEM[dim],DqEM[dim];

     val_type da,db;
     val_type *Gm,GmSun,Gma[nbody];

     val_type GQN[dim];
     val_type sub1[dim],sub2[dim],sub3[dim];

     Pkepler_sys *Pkepler;
     val_type *mu,*k2;

/* ----------- implementation  ---------------------------------------*/

//     printf("DRR2B ebaluazioa\n");

     nd=neq/2;
     Gm=params->rpar;
     Pkepler=&params->Pkepler;
     mu=Pkepler->Mu;
     k2=Pkepler->K2;

     iM=(nbody-1);  	       // Moon
     iM1=iM*dim;
     iM2=nd+iM1;
     iE=(nbody-2);            // Earth
     iE1=iE*dim;
     iE2=nd+iE1;

     GmSun=Gm[0];
     for (i=0; i<nbody; i++) {Gma[i]=Gm[i+1];}
     for (id=0; id<dim; id++) GQN[id]=-(u[iM1+id]/Gma[iE])*Gma[iM];


     for (id=0; id<dim; id++)
     {
          DqEM[id]=(u[iM1+id]/Gma[iE])*(Gma[iM]+Gma[iE]);
          qEM[id]=-(u[iE1+id]+u[iM1+id]);

     }


     for (i=0; i<nbody; i++)
     {
          i1=i*dim;
          i2=nd+i1;
          for (id=0; id<dim; id++)
          {
               dR[i1+id]=0.;
               dR[i2+id]=0.;
          }
     }

     /* First DQi  */


     for (i=0; i<(nbody-1); i++)
     {
          i1=i*dim;
          i2=nd+i1;

          for (j=0; j<(nbody-1); j++)
          {
               j1=j*dim;
               j2=nd+j1;
  
               if (j!=i)
                   for (id=0; id<dim; id++)
                        dR[i1+id]+=(u[j2+id]*mu[j]);
           }

          for (id=0; id<dim; id++)
                 dR[i1+id]*=1/GmSun;

     }


     /* Second DVi  */


     for (i=0; i<(nbody-2); i++)
     {
           i1=i*dim;
           i2=nd+i1;

           for (j=0; j<(nbody-2); j++)
           {
                j1=j*dim;
                j2=nd+j1;

                if (j!=i)
                {
                    da=0;
                    for (id=0; id<dim; id++)
                    {
                         qa[id]=(u[i1+id]-u[j1+id]);
                         da+=qa[id]*qa[id];
                    }

                    da=SQRT(da)*da;

                    for (id=0; id<dim; id++)
                         dR[i2+id]-=Gma[j]*qa[id]/da;
                }
           }


           da=0.;
           db=0.;
           for (id=0; id<dim; id++)
           {
                 qa[id]=u[i1+id]-u[iE1+id]-u[iM1+id];
                 qb[id]=u[i1+id]-u[iE1+id]-GQN[id];     
                 da+=qa[id]*qa[id];
                 db+=qb[id]*qb[id];
           }

           da=SQRT(da)*da;
           db=SQRT(db)*db;

           for (id=0; id<dim; id++)
           {
                 aux1=Gma[iM]/da*qa[id]+Gma[iE]/db*qb[id];
                 aux2=Gma[iM]*(qa[id]/da-qb[id]/db);
                 dR[i2+id]-=aux1;
                 dR[iE2+id]+=Gma[i]*aux1;
                 dR[iM2+id]+=Gma[i]*aux2;
           }

           for (id=0; id<dim; id++)
                 dR[i2+id]=(dR[i2+id]/mu[i])*Gma[i];
     }

     /* Dvi i=N-1, i=N */

     mysubstract(&u[iE1],&u[iM1],sub1);
     mysubstract(&u[iE1],GQN,sub2);

     mysubstract(qEM,DqEM,sub3);


     for (id=0; id<dim; id++)
     {
          dR[iE2+id]-=GmSun*(Gma[iE]*sub2[id]+Gma[iM]*sub1[id]);
          dR[iM2+id]+=-GmSun*(Gma[iM]*sub3[id]);
          dR[iE2+id]*=1/mu[iE];
          dR[iM2+id]*=1/mu[iM];

     }


     for (i=0; i<nbody; i++)
     {
          i1=i*dim;
          i2=nd+i1;
   
          da=0.;
          for (id=0; id<dim; id++) da+=u[i1+id]*u[i1+id];
          da=SQRT(da)*da;
                 
          for (id=0; id<dim; id++)
               dR[i2+id]-=k2[i]*u[i1+id]/da;

     }


#ifdef MDEBUG

    double myaux;
    printf("\nDRR2B\n");
    printf("q :");    
    for (i=0; i<neq/2;i++)
    {
         myaux=dR[i];
         printf("%lg,", myaux);
    }
    printf("\n");
    printf("v :");
    for (i=neq/2; i<neq;i++) 
    {
        myaux=dR[i]; 
        printf("%lg,", myaux);
    }
    printf("\n");


#endif

     return ;

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



