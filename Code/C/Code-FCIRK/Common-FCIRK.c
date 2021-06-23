/*----------------------------------------------------------------------------*/
/*									      */
/*    Common-FCIRK.c							      */
/*									      */
/*	Functions: 							      */
/*	 InitStat()							      */
/*	 NormalizedDistance()						      */
/*	 IRKstep_fixed()						      */
/*	 FP_Iteration()		         			      */
/*	 Summation()							      */
/*       deltafun	                                                    */
/*       Ordinary_stepQ                                                    */
/*       Num_steps                                                         */
/*	 Main_FCIRK()							      */
/*       								      */
/* ---------------------------------------------------------------------------*/


#include <Common-FCIRK.h>
#include <float.h>

FILE *myfile2;


/******************************************************************************/
/* 					   				      */
/* InitStat 		         	  				      */
/*                                        				      */
/*									      */
/******************************************************************************/

void InitStat (ode_sys *system, tcoeffs *gsmethod,
               tcache_stat *cache_stat, tcache_vars *cache_vars )
{

     int i,is,in,neq,ns;

     ns=gsmethod->ns;
     neq=system->neq;

     cache_stat->laststep = false;
     cache_stat->stepcount = 0;
     cache_stat->hstepcount = 0;
     cache_stat->itcount=0;
     cache_stat->totitcount=0;
     cache_stat->totitcountzero=0;
     cache_stat->maxitcount=0;
     cache_stat->itzero=0;
     cache_stat->fcn=0;
     cache_stat->convergence=SUCCESS;
     cache_stat->nout=0;
     cache_stat->MaxDE=0.;

     cache_vars->z = (val_type *)malloc(neq*(ns)*sizeof(val_type));
     cache_vars->li = (val_type *)malloc(neq*ns*sizeof(val_type));
     cache_vars->fz = (val_type *)malloc(neq*ns*sizeof(val_type));
     cache_vars->zold = (val_type *)malloc(neq*ns*sizeof(val_type));
     cache_vars->DMin = (val_type *)malloc(neq*ns*sizeof(val_type));
     cache_vars->delta = (val_type *)malloc(neq*sizeof(val_type));
     cache_vars->avg_delta = (val_type *)malloc(neq*sizeof(val_type));

     for (is=0; is<ns; is++)
     {
          in=is*neq;
          for (i=0; i<neq; i++)
          {
               cache_vars->li[in+i]=0.;
               cache_vars->fz[in+i]=0.;
          }
     }


    return;

}


/******************************************************************************/
/* 					   				      */
/* InitStat_high	         	  				      */
/*                                        				      */
/*									      */
/******************************************************************************/

void InitStat_high (ode_sys_high *system, tcoeffs_h *coeffs,
                    tcache_vars_high *cache_vars )
{

     int i,is,in,neq,ns;

     ns=coeffs->ns;
     neq=system->neq;

     cache_vars->z = (highprec *)malloc(neq*(ns)*sizeof(highprec));
     cache_vars->li = (highprec *)malloc(neq*ns*sizeof(highprec));
     cache_vars->fz = (highprec *)malloc(neq*ns*sizeof(highprec));
     cache_vars->zold = (highprec *)malloc(neq*ns*sizeof(highprec));
     cache_vars->DMin = (highprec *)malloc(neq*ns*sizeof(highprec));


     for (is=0; is<ns; is++)
     {
          in=is*neq;
          for (i=0; i<neq; i++)
          {
               cache_vars->li[in+i]=0.;
               cache_vars->fz[in+i]=0.;
          }
     }


    return;

}




/******************************************************************************/
/* 					   				      */
/* NormalizedDistance 	         	  				      */
/*                                        				      */
/*									      */
/******************************************************************************/

val_type NormalizedDistance ( int neq, int ns,
                              toptions *options,  val_type *z,
                              val_type *zold)
{

#define BASE val_type
#define HIGH 0
#include <source_NormalizedDistance.h>
#undef BASE
#undef HIGH

}


val_type NormalizedDistance_high ( int neq, int ns,
                                   toptions *options,  highprec *z,
                                   highprec *zold)
{

#define BASE highprec
#define HIGH 1
#include <source_NormalizedDistance.h>
#undef BASE
#undef HIGH

}

/******************************************************************************/
/* 									      */
/*   IRKstep_fixed							      */
/* 									      */
/******************************************************************************/


void IRKstep_fixed         (  ode_sys *system,  solution *u,
                              val_type tn, int ii, val_type h,
                              toptions *options,  tcoeffs *method,
                              tcache_stat *cache_stat, tcache_vars *cache_vars)
{

#define BASE val_type
#define HIGH 0
#include <source_IRKstep_fixed.h>
#undef BASE
#undef HIGH

}



void IRKstep_fixed_high   (   ode_sys_high *system,  solution *u,
                              highprec tn, int ii, highprec h,
                              toptions *options,  tcoeffs_h *method,
                              tcache_stat *cache_stat, tcache_vars_high *cache_vars)
{

#define BASE highprec
#define HIGH 1


#include <source_IRKstep_fixed.h>
#undef BASE
#undef HIGH

}


/******************************************************************************/
/* 									      */
/*   IRKstep_adaptive							      */
/* 									      */
/******************************************************************************/


void IRKstep_adaptive  (tode_sys *ode_system,solution *u, 
                        val_type tn, int ii, val_type h,toptions *options,
                        tmethod *method, 
                        tcache_stat *cache_stat, 
                        tcache_vars *cache_vars, tcache_vars_high *cache_vars_high)

{

#ifdef PARALLEL
     int extern thread_count;
#endif

/* ----- First initializations -----------------------------------------------*/
     int neq,ns;

     neq=ode_system->system.neq;
     ns=method->coeffs.ns;

/*------ Declarations --------------------------------------------------------*/

     int i,is,isn,in,js,jsn,ism,ismi;
     int k,ki;
     bool iter0;
     int D0;
     double difftest;
     val_type *z,*zold,*fz,*li,*DMin,*coef;
     val_type *ttau;    
     val_type sum; 
     parameters *params;
     Pkepler_sys *Pkepler;

     highprec tj,hj;

/* ----- Implementation ------------------------------------------------------*/


     z=cache_vars->z;
     zold=cache_vars->zold;
     fz=cache_vars->fz;
     li=cache_vars->li;
     DMin=cache_vars->DMin;
     ttau=method->coeffs.ttau;
     coef=method->coeffs.nu;

     params=&ode_system->system.params;
     Pkepler=&params->Pkepler;
       

/* --- Interpolation ---------------------------------------------------------*/


     if (cache_stat->interpolate !=true)


        for (is = 0; is<ns; is++)
        {
            in=is*neq;
            for (i = 0; i<neq; i++) z[in+i]=u->uul[i];
        }

     else
     {
       for (is = 0; is<ns; is++)
       {

          in=is*neq;
          KeplerFlowAll2 (neq,Pkepler->keplerkop,&z[in],h,params);

          for (i = 0; i<neq; i++)
          {
              zold[in+i]=z[in+i];
          }
       }

       for (is = 0; is<ns; is++)
       {
          in=is*neq;
          for (i = 0; i<neq; i++)
          {
              sum=0.;
              for (js = 0; js<ns; js++)
              {
                   jsn=(ns+1)*is+js;
                   sum+=zold[neq*js+i]*(coef[jsn]);
              }

             jsn=(ns+1)*is+ns;
             z[in+i]=sum+u->uul[i]*coef[jsn];

          }
       }
     } 


/*---- Fixed point iteration ------------------------------------------------*/


     iter0=true;
     D0=0;
     for (is=0; is<ns; is++)
     {
          isn=neq*is;
          for (i=0; i<neq; i++)
          {
               in=isn+i;
               DMin[in]=INF;
          }
     }

     cache_stat->itcount=0;
     ism=ii*ns;

#ifdef PARALLEL
#      pragma omp parallel for num_threads(thread_count) private(isn,ismi)
#endif
     for (is = 0; is<ns; is++)
     {
           isn=neq*is;
           ismi=ism+is;
           cache_stat->fcn++;         
           ode_system->system.f(neq,tn+method->coeffs.hc[is],ttau[ismi],&z[isn],&fz[isn],params);
           for (i=0; i<neq; i++) li[isn+i]=fz[isn+i]*method->coeffs.hb[is];
     }

     cache_stat->itcount++;

     while (iter0 && cache_stat->itcount<MAXIT)
     {
           FP_Iteration (&ode_system->system,u,tn,ii,h,&method->coeffs,cache_stat,cache_vars,&D0, &iter0);      
           cache_stat->itcount++;

     }


     if (cache_stat->itcount==MAXIT)
     {
           printf("Break: step(MAXIT)=%i\n",cache_stat->itcount);
           cache_stat->convergence=FAIL;
     }


     if (D0<(ns*neq))
     {
           difftest=NormalizedDistance(neq,ns,options,z,zold);
           if (difftest>1.)
           {
                cache_stat->convergence=FAIL;
                printf("Lack of convegence of Fixed point iteration:\
                        step=%i,iteration=%i,",
                        cache_stat->stepcount,cache_stat->itcount);
                printf("difftest=%lg\n",difftest);
           }
     }
     else
     {
          (cache_stat->totitcountzero)+=(cache_stat->itcount);
          cache_stat->itzero++;
   
     }


     if (cache_stat->convergence!=FAIL)
     {
           deltafun(&method->coeffs,&ode_system->system,cache_vars,&cache_vars->delta[0]);
           k=Num_steps(&method->coeffs,&ode_system->system,cache_stat,cache_vars);
           if (Ordinary_stepQ(&ode_system->system,h,k,cache_stat,cache_vars)==true)
           {
             Summation (&method->coeffs,u,&ode_system->system,options,cache_vars);
           }
           else
           { 
             cache_stat->hstepcount++;
             hj=h/k;
             AssignGsCoefficients_high(DIR_COEFF,&method->coeffs_h,hj,k);
             UpdateGsCoefficients_high(&method->coeffs_h,hj,k);
             cache_stat->interpolate=false;
             for (ki=0; ki<k; ki++)
             {
                  tj=tn+ki*hj;   
                  IRKstep_fixed_high (&ode_system->system_h,u,tj,ki,hj,options,&method->coeffs_h,cache_stat,cache_vars_high);                
             }
           }    
     }

     return;

}




/******************************************************************************/
/*									      */
/*      FP_Iteration: Fixed Point Iteration	         	      	      */
/*									      */
/******************************************************************************/


int FP_Iteration     ( ode_sys *system,  solution *u,  val_type tn,
                         int ii,val_type h,  tcoeffs *method,
                         tcache_stat *cache_stat, tcache_vars *cache_vars,
                         int *D0, bool *iter0)
{
#define BASE val_type
#define KFABS(x) FABS(x)
#define HIGH 0
#include <source_FP_iteration.h>
#undef BASE
#undef KFABS
#undef HIGH

}


int FP_Iteration_high  ( ode_sys_high *system,  solution *u,  highprec tn,
                         int ii,highprec h,  tcoeffs_h *method,
                         tcache_stat *cache_stat, tcache_vars_high *cache_vars,
                         int *D0, bool *iter0)
{
#define BASE highprec
#define KFABS(x) FABS_high(x)
#define HIGH 1
#include <source_FP_iteration.h>
#undef BASE
#undef KFABS
#undef HIGH

}

/******************************************************************************/
/* 									      */
/*    Summation								      */
/* 									      */
/* 									      */
/******************************************************************************/

void Summation                ( tcoeffs *gsmethod,
                                solution *u, ode_sys *system,
                                toptions *options, tcache_vars *cache_vars)

{

#define BASE val_type
#define HIGH 0
#include <source_Summation.h>
#undef BASE
#undef HIGH

}


void Summation_high            ( tcoeffs_h *gsmethod,
                                solution *u, ode_sys_high *system,
                                toptions *options, tcache_vars_high *cache_vars)

{

#define BASE highprec
#define HIGH 1
#include <source_Summation.h>
#undef BASE
#undef HIGH

}


/******************************************************************************/
/* 									      */
/*   deltafun    							      */
/* 									      */
/* 									      */
/******************************************************************************/

void deltafun     ( tcoeffs *gsmethod, ode_sys *system,
                    tcache_vars *cache_vars, val_type *delta)

{


/* ----- First initializations -----------------------------------------------*/

     int neq,ns;

     neq=system->neq;
     ns=gsmethod->ns;

/*------ Declarations --------------------------------------------------------*/

     int i,is,isn;
     val_type *fz,*d;
     fz=cache_vars->fz;
     d=gsmethod->d;

/* ----- Implementation ------------------------------------------------------*/


    for (i = 0; i<neq; i++)
    {
        delta[i]=fz[i]*d[0];

        for (is=1; is<ns; is++)
        {
                isn=neq*is;
                delta[i]+=fz[isn+i]*d[is];
        }
    }


#ifdef DELDEBUG

   double aux;

   aux=0.;
   for (i=0; i<neq; i++) aux+=delta[i]*delta[i];
   printf("Norm(delta)=%lg\n",sqrt(aux));
 

#endif


    return;

}



/******************************************************************************/
/* 									      */
/*    Ordinary_stepQ							      */
/* 									      */
/* 									      */
/******************************************************************************/

bool Ordinary_stepQ            ( ode_sys *system, val_type h, int k,
                                 tcache_stat *cache_stat, tcache_vars *cache_vars)

{


/* ----- First initializations -----------------------------------------------*/

     int neq;
     neq=system->neq;
 

/*------ Declarations --------------------------------------------------------*/

     int i,m;
     bool res;

     val_type *delta,*avg_delta;
     delta=cache_vars->delta;
     avg_delta=cache_vars->avg_delta;

/* ----- Implementation ------------------------------------------------------*/

    res=true;
    m=cache_stat->stepcount+1;

    if (m>1)
    {
       for (i=0; i<neq; i++)
       {
           if (FABS(delta[i])>NU*avg_delta[i]) 
           {
             res=false;
#ifdef IOUT
             double lag;
             lag=cache_stat->stepcount*h;
             fprintf(myfile2,"Non-Ordinary step!!. step=%i, t=%lg, gorputza=%i,k=%i\n",cache_stat->stepcount,lag,i,k);
#endif
           }
           avg_delta[i]=(avg_delta[i]*(m-1)+FABS(delta[i]))/m;           
       }
    }
    else
       for (i = 0; i<neq; i++) avg_delta[i]=FABS(delta[i]);

    

#ifdef MDEBUG

   double aux;

   aux=0.;
   for (i=0; i<neq; i++) aux+=avg_delta[i]*avg_delta[i];
   printf("Norm(avg_delta)=%lg\n",sqrt(aux));

 
#endif



    return(res);

}




/******************************************************************************/
/* 									      */
/*    Num_steps								      */
/* 									      */
/* 									      */
/******************************************************************************/

int Num_steps     ( tcoeffs *gsmethod,
                    ode_sys *system,
                    tcache_stat *cache_stat, 
                    tcache_vars *cache_vars)
{


/* ----- First initializations -----------------------------------------------*/

     int neq,ns;

     neq=system->neq;
     ns=gsmethod->ns;
 
/*------ Declarations --------------------------------------------------------*/

     int i,k,m;
     val_type max,p;
     val_type *delta,*avg_delta;
     delta=cache_vars->delta;
     avg_delta=cache_vars->avg_delta;

/* ----- Implementation ------------------------------------------------------*/

     m=cache_stat->stepcount+1;

    if (m>1)
    {
       max=0.;
       p=1./(ns-1);
       for (i = 0; i<neq; i++)
          max=FMAX(max,pow(FABS(delta[i])/avg_delta[i],p));

       k=ceil(max);
    }
    else
      k=1;

   return(k);


#ifdef ESTDEBUG

   double a;

   a=max;
   printf("Max=%lg,k=%i\n",a,k);

 
#endif


}






/******************************************************************************/
/* 									      */
/*   Main-FCIRK	  							      */
/*   									      */
/* 									      */
/******************************************************************************/

void Main_FCIRK
(val_type t0, val_type t1, val_type h,
 tmethod *method, solution *u, tode_sys *ode_system, toptions *options,
 tcache *cache)

{

/* ----- First initializations -----------------------------------------------*/

     int neq;

     parameters *params;

     neq=ode_system->system.neq;
     params=&ode_system->system.params;


/*------ Declarations --------------------------------------------------------*/

     FILE *myfile;

     int istep,nstep;
     val_type tn;
     val_type *fz;

     fz=cache->cache_vars.fz;

/* ----- Implementation ------------------------------------------------------*/

#ifdef IOUT
     myfile = fopen(options->filename,"wb");
     myfile2 = fopen("./Data/NonOrdinarySteps.txt","w");
#endif

     tn=t0;
     nstep=round((t1-t0)/h);
     cache->cache_stat.execution=SUCCESS;

#ifdef IOUT
     options->TheOutput(&ode_system->system,&method->coeffs,tn,h,u,&cache->cache_stat,params,options,myfile);
#endif

     ode_system->system.StartFun(neq,t0,h,u,params);
     cache->cache_stat.interpolate=false;

     for(istep=0; (istep<nstep); istep++)
     {
          if (options->adaptive==true)
          {
             IRKstep_adaptive (ode_system,u,tn,0,h,options,
                               method,
                               &cache->cache_stat,&cache->cache_vars,&cache->cache_vars_h);
          }
          else
          {  
             IRKstep_fixed (&ode_system->system,u,tn,0,h,options,&method->coeffs,&cache->cache_stat,&cache->cache_vars);     

          }
     
          if (h<7) cache->cache_stat.interpolate=true;

          if (cache->cache_stat.convergence==FAIL || isnan(fz[0]))
          {
               cache->cache_stat.execution=FAIL;
               nstep=istep;

          }
          else
          {
               tn=t0+(istep+1)*h;

               cache->cache_stat.stepcount++;
               (cache->cache_stat.totitcount)+=(cache->cache_stat.itcount);
               if ((cache->cache_stat.itcount)>(cache->cache_stat.maxitcount))
                  (cache->cache_stat.maxitcount)=(cache->cache_stat.itcount);

               ode_system->system.ProjFun(neq, tn, h, u, params);
#ifdef IOUT
               options->TheOutput(&ode_system->system,&method->coeffs,tn,h,u,&cache->cache_stat,params,options,myfile);
#endif
               if (tn+h>=t1 && cache->cache_stat.laststep==false)
                    cache->cache_stat.laststep=true;

          }
     }

#ifdef DEBUG
     if (cache->cache_stat.execution==SUCCESS) printf("Success execution!!\n");
     else printf("Execution Fail!!\n");
#endif

#ifdef IOUT
     fclose(myfile);
     fclose(myfile2);
#endif

     return;

}




