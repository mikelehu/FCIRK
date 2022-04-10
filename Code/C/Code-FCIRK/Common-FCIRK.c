/*----------------------------------------------------------------------------*/
/*									      */
/*    Common-FCIRK.c							      */
/*									      */
/*	Functions: 							      */
/*	 InitStat()							      */
/*       add2 and sub2				       		      */
/*	 NormalizedDistance()						      */
/*	 IRKstep_fixed()						      */
/*	 FP_Iteration()		         			      */
/*	 Summation()							      */
/*	 UnSummation()							      */
/*       deltafun	                                                    */
/*       Ordinary_stepQ                                                    */
/*       Num_steps                                                         */
/*       Rad_Convergence                                                    */
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
     int dim,nbody;

     ns=gsmethod->ns;
     neq=system->neq;
     
     dim=3;
     nbody=neq/(2*dim);

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
     cache_vars->ux1 = (val_type *)malloc(neq*sizeof(val_type));
     cache_vars->ux2 = (val_type *)malloc(neq*sizeof(val_type));
     cache_vars->uB = (val_type *)malloc((neq+2*dim)*sizeof(val_type));
     cache_vars->KK = (val_type *)malloc((nbody+1)*sizeof(val_type));
     
     cache_vars->RC=0.;
     cache_vars->MaxRC=0.;
     cache_vars->MinRC=1.e6;
     cache_vars->RCSum=0.;
     cache_vars->RCSum2=0.;

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


/*************************************************************************************/
/* 					   				                */
/* add2  (Dekker)              	  			                       */
/*                                        					       */
/*                                        					       */
/*     add2 calculates the doublelength sum of (x, xx) and (y, yy),                  */
/*    the result being (z, zz) ;                                       	        */
/*									               */
/************************************************************************************/
void add2 (val_type x, val_type xx, val_type y, val_type yy,
           val_type *z, val_type *zz)
{   
/*------ declarations ---------------------------------------------------*/

     val_type  s,r;

/*------ implementaion --------------------------------------------------*/   
 
     r=x+y;
     if (FABS(x) > FABS(y))
        s=x-r+y+yy+xx;
     else
        s=y-r+x+xx+yy;
        
     *z=r+s;
     *zz=r-*z+s;                             
 
     return;

 }


void sub2 (val_type x, val_type xx, val_type y, val_type yy,
           val_type *z, val_type *zz)
{   
/*------ declarations ---------------------------------------------------*/

     val_type  s,r;

/*------ implementaion --------------------------------------------------*/   
 
     r=x-y;
     if (FABS(x) > FABS(y))
        s=x-r-y-yy+xx;
     else
        s=-y-r+x+xx-yy;
        
     *z=r+s;
     *zz=r-*z+s;                             
 
     return;

 }


/******************************************************************************/
/* 					   				        */
/* Rmdigits      (31-10-2021)        	  				        */
/*                                        				        */
/*									        */
/******************************************************************************/

val_type Rmdigits (val_type x,val_type r)
{
  
#define BASE val_type
#define HIGH 0
#include <source_Rmdigits.h>
#undef BASE
#undef HIGH

}


/******************************************************************************/
/* 					   				        */
/* Rmdigits_high      (31-10-2021)  	  				        */
/*                                        				        */
/*									        */
/******************************************************************************/

highprec Rmdigits_high (highprec x,highprec r)
{

#define BASE highprec
#define HIGH 1
#include <source_Rmdigits.h>
#undef BASE
#undef HIGH
   


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
     cache_vars->ux1 = (highprec *)malloc(neq*(ns)*sizeof(highprec));
     cache_vars->ux2 = (highprec *)malloc(neq*(ns)*sizeof(highprec));


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
     
     val_type R,dR;             //2021-12-10  
       

/* --- Interpolation ---------------------------------------------------------*/


     if (cache_stat->interpolate !=true)


        for (is = 0; is<ns; is++)
        {
            in=is*neq;
            for (i = 0; i<neq; i++) z[in+i]=u->uu[i];
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
             z[in+i]=sum+u->uu[i]*coef[jsn];

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
           ode_system->system.f(neq,tn+method->coeffs.hc[is],ttau[ismi],
                                &z[isn],&fz[isn],params,cache_vars);
           for (i=0; i<neq; i++) li[isn+i]=fz[isn+i]*method->coeffs.hb[is];
     }

     cache_stat->itcount++;

     while (iter0 && cache_stat->itcount<MAXIT)
     {
           FP_Iteration (&ode_system->system,u,tn,ii,h,
                         options,&method->coeffs,cache_stat,cache_vars,&D0, &iter0);      
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

           Summation (&method->coeffs,u,&ode_system->system,options,cache_vars);
          
           cache_vars->RC=ode_system->system.rc_fun(&ode_system->system,cache_vars,u); 
           cache_vars->RCSum+=cache_vars->RC;
           cache_vars->RCSum2+=(cache_vars->RC)*(cache_vars->RC);
           if (cache_vars->RC>cache_vars->MaxRC) cache_vars->MaxRC=cache_vars->RC;
           if (cache_vars->RC<cache_vars->MinRC) cache_vars->MinRC=cache_vars->RC;
           
           R=cache_vars->RCSum/(cache_stat->stepcount+1);
           dR=SQRT(cache_vars->RCSum2/(cache_stat->stepcount+1)-R*R);

           
           if (NU*cache_vars->RC<(R-dR))
           {
              UnSummation (&method->coeffs,u,&ode_system->system,options,cache_vars);
              k=CEIL(R/cache_vars->RC);              
              cache_stat->hstepcount++;
              hj=h/k;
              AssignGsCoefficients_high(DIR_COEFF,&method->coeffs_h,hj,k);
              UpdateGsCoefficients_high(&method->coeffs_h,hj,k);
              cache_stat->interpolate=false;
#ifdef IOUT
              double lag;
              lag=tn;
              fprintf(myfile2,"Non-Ordinary step!!. step=%i, t=%lg,k=%i\n",cache_stat->stepcount,lag,k);
#endif
              for (ki=0; ki<k; ki++)
              {
                  tj=tn+ki*hj;   
                  IRKstep_fixed_high (&ode_system->system_h,u,tj,ki,hj,options,
                                      &method->coeffs_h,cache_stat,cache_vars_high);                
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
                         int ii,val_type h, toptions *options,  tcoeffs *method,
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
                         int ii,highprec h, toptions *options, tcoeffs_h *method,
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
/*    Summation							      */
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
#define LOW val_type
#define BASE highprec
#define HIGH 1
#include <source_Summation.h>
#undef LOW
#undef BASE
#undef HIGH

}



/******************************************************************************/
/* 									      */
/*    UnSummation							      */
/* 									      */
/* 									      */
/******************************************************************************/

void UnSummation                ( tcoeffs *gsmethod,
                                solution *u, ode_sys *system,
                                toptions *options, tcache_vars *cache_vars)

{

#define BASE val_type
#define HIGH 0
#include <source_UnSummation.h>
#undef BASE
#undef HIGH

}


void UnSummation_high            ( tcoeffs_h *gsmethod,
                                solution *u, ode_sys_high *system,
                                toptions *options, tcache_vars_high *cache_vars)
{
#define LOW val_type
#define BASE highprec
#define HIGH 1
#include <source_UnSummation.h>
#undef LOW
#undef BASE
#undef HIGH

}



/******************************************************************************/
/* 									        */
/*   Rad_Convergence  							        */
/* 	Lower bound of the radius of convergence of the N-body                */
/* 									        */
/******************************************************************************/


val_type Rad_Convergence1 (ode_sys *system,
                          tcache_vars *cache_vars, solution *U)
{
/*  15-body problem */

  return Rad_Convergence (system,cache_vars,U,ChangeHeltoBar_EMB,false);

}

val_type Rad_Convergence2 (ode_sys *system,
                          tcache_vars *cache_vars, solution *U)
{
/*  16-body problem : Moon is incluided*/

  return Rad_Convergence (system,cache_vars,U,ChangeHeltoBar_Moon,false);

}


val_type Rad_Convergence3 (ode_sys *system,
                          tcache_vars *cache_vars, solution *U)
{
/*  16-body problem : Moon is excluided*/

  return Rad_Convergence (system,cache_vars,U,ChangeHeltoBar_Moon,true);

}


/******************************************************************************/

val_type Rad_Convergence (ode_sys *system,
                          tcache_vars *cache_vars, solution *U,
                          void ChangeHeltoBar(),
                          bool exc_Moon)

{

/*------ Declarations --------------------------------------------------------*/

     int dim,neqH,nbodyH;
     int neqB,nbodyB,ndB;
     int i,i1,i2,j,j1,j2;
     int iE,iM;
     
     val_type xi,yi,zi,vxi,vyi,vzi;
     val_type xij,yij,zij,vxij,vyij,vzij;
     val_type norm2qij,normqij,normvij;
     val_type a,b,L,Lij;
     
     parameters *params;
     Pkepler_sys *Pkepler;
     val_type *KK,*u,*Gm,*mu;               

/* ----- Implementation ------------------------------------------------------*/

     dim=3;
     neqH=system->neq;
     nbodyH=neqH/(2*dim);
     
     neqB=neqH+2*dim;
     nbodyB=neqB/(2*dim);
     ndB=neqB/2;
     
     iE=nbodyB-2;  // Earth baryc. index n-1
     iM=nbodyB-1;  // Moon baryc,  index n
     
     KK=cache_vars->KK;         
     u=cache_vars->uB;
     params=&system->params;
     Pkepler=&params->Pkepler;
     Gm=params->rpar;
     mu=Pkepler->Mu;

     
     ChangeHeltoBar (neqH, nbodyH, u, U, Gm, mu);
        
     
     for (i=0; i<nbodyB; i++) KK[i]=0.;
     
     for (i=0; i<nbodyB; i++)
     {
        i1=i*dim;
        xi=u[i1];
        yi=u[i1+1];
        zi=u[i1+2];
        
        for (j=i+1;j<nbodyB;j++)
        {
            j1=j*dim;
            xij=xi-u[j1];
            yij=yi-u[j1+1];
            zij=zi-u[j1+2];
            norm2qij=xij*xij+yij*yij+zij*zij;
            KK[i]+=Gm[j]/norm2qij; 
            KK[j]+=Gm[i]/norm2qij;

        }
     }
     
     L=0.;
     
     for (i=0; i<nbodyB; i++)
     {
        i1=i*dim;
        i2=ndB+i1;
        xi=u[i1];
        yi=u[i1+1];
        zi=u[i1+2];
        vxi=u[i2];
        vyi=u[i2+1];
        vzi=u[i2+2];     

        for (j=i+1;j<nbodyB;j++)
        {
            j1=j*dim;
            j2=ndB+j1;
            xij=xi-u[j1];
            yij=yi-u[j1+1];
            zij=zi-u[j1+2];
            vxij=vxi-u[j2];
            vyij=vyi-u[j2+1];
            vzij=vzi-u[j2+2];
            norm2qij=xij*xij+yij*yij+zij*zij;
            normqij=SQRT(norm2qij);
            normvij=SQRT(vxij*vxij+vyij*vyij+vzij*vzij);
	    a=normvij/normqij;
	    b=(KK[i]+KK[j])/normqij;
	    if (!(exc_Moon==true && i==iE && j==iM))
	    {
	        Lij=a+SQRT(a*a+4./7*b);
	        L=FMAX(L,Lij);
	    }
        }
     }   


#ifdef XDEBUG

   double aux;

   aux=2./(7*L);
   printf("Radius of convergence%lg\n",aux);
 
#endif

  
  return (2./(7*L));

}






/******************************************************************************/
/* 									      */
/*   Main-FCIRK	  						      */
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
     
     solution uj2;
     highprec errm,errj[neq];
     val_type zz[2];
     val_type h2,tj2;
     
     tcoeffs coeffs2; 
     tcache_stat cache_stat2;
     tcache_vars cache_vars2;

     Pkepler_sys *Pkepler;
     Pkepler=&params->Pkepler;

/*------ Declarations --------------------------------------------------------*/

     FILE *myfile;
     FILE *myfileRC,*myfileER;
//     *myfileEL ,*myfileERZ;

     int istep,nstep;
     int i,ii,l;
     val_type tn;
     val_type *fz;

     fz=cache->cache_vars.fz;

/* ----- Implementation ------------------------------------------------------*/

#ifdef IOUT
     myfile = fopen(options->filename,"wb");
     myfile2 = fopen("./Data/NonOrdinarySteps.txt","w");
     myfileRC = fopen("./Data/RConvergence.bin","wb");
     myfileER = fopen("./Data/LocalError.bin","wb");
     
     fprintf(myfile2,"Non-Ordinary step!!\n");
#endif


     tn=t0;
     nstep=round((t1-t0)/h);
     cache->cache_stat.execution=SUCCESS;
     
     if (options->errorsL==true)
        for (i=0; i<neq; i++) errj[i]=0.;
        
     uj2.uu = (val_type *)malloc(neq*sizeof(val_type));
     uj2.ee = (val_type *)malloc(neq*sizeof(val_type));        
     
     coeffs2.ns=8;
     h2=h/2.; 
     MallocGsCoefficients (DIR_COEFF,&coeffs2,KMAX);
     AssignGsCoefficients(DIR_COEFF,&coeffs2,h2,1);
     
     InitStat(&ode_system->system,&coeffs2,&cache_stat2,&cache_vars2);
     cache_stat2.interpolate=false;   
     
     for (i=0; i<neq; i++)
     {
       uj2.uu[i]=u->uu[i];
       uj2.ee[i]=u->ee[i];
      }   
          
#ifdef IOUT
     options->TheOutput(&ode_system->system,&method->coeffs,tn,h,u,
                        &cache->cache_stat,&cache->cache_vars,params,options,
                        myfile,myfileRC,myfileER,errj);
#endif

     ode_system->system.StartFun(neq,t0,h,u,params);
     cache->cache_stat.interpolate=false;          
     
     
     solution ux1;
     ux1.uu = (val_type *)malloc((neq)*sizeof(val_type)); 
     ux1.ee = (val_type *)malloc((neq)*sizeof(val_type)); 

     for(istep=0; (istep<nstep); istep++)
     {
     
          tj2=tn;
          

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
          
               if (options->errorsL==true)
               {        
                
                  for (ii=0; ii<2; ii++)
                  {

			KeplerFlowAll_high (neq,Pkepler->keplerkop,&uj2, h2/2.,params); 		                      
                       IRKstep_fixed (&ode_system->system,&uj2,tj2,0,h2,options,&coeffs2,&cache_stat2,&cache_vars2); 
			KeplerFlowAll_high (neq,Pkepler->keplerkop,&uj2, h2/2.,params); 

                       tj2+=h2;
                  }
                  
                 for(i=0; i<neq; i++)
                 {
                   ux1.uu[i]=u->uu[i];
                   ux1.ee[i]=u->ee[i];
                 }
                  
                 KeplerFlowAll_high (neq,Pkepler->keplerkop,&ux1, h/2.,params); 
                                         
                               
                  for (l=0; l<neq; l++)
                  {
                       sub2(ux1.uu[l],ux1.ee[l],uj2.uu[l],uj2.ee[l],&zz[0],&zz[1]);
                       errm=zz[0];
                       errm+=zz[1];
                       if (errj[l]<FABS_high(errm)) errj[l]=FABS_high(errm); 
                  }
               
             }          
             
             
               tn=t0+(istep+1)*h;

               cache->cache_stat.stepcount++;
               (cache->cache_stat.totitcount)+=(cache->cache_stat.itcount);
               if ((cache->cache_stat.itcount)>(cache->cache_stat.maxitcount))
                  (cache->cache_stat.maxitcount)=(cache->cache_stat.itcount);
                  
               for (i=0; i<neq; i++)
               {
                  uj2.uu[i]=ux1.uu[i];
                  uj2.ee[i]=ux1.ee[i];
               }  

               ode_system->system.ProjFun(neq, tn, h, u, params);
#ifdef IOUT
               options->TheOutput(&ode_system->system,&method->coeffs,tn,h,u,
                                  &cache->cache_stat,&cache->cache_vars,params,options,
                                  myfile,myfileRC,myfileER,errj);
#endif
               if (FABS(tn+h)>=FABS(t1) && cache->cache_stat.laststep==false)
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
     fclose(myfileRC);
     fclose(myfileER);
#endif

     free(uj2.uu);
     free(uj2.ee);
     
     free(ux1.uu);   
     free(ux1.ee);
         
     free(coeffs2.m);
     free(coeffs2.a);
     free(coeffs2.b);
     free(coeffs2.hb);
     free(coeffs2.c);
     free(coeffs2.hc);
     free(coeffs2.nu);
     free(coeffs2.ttau);

     free(cache_vars2.z);
     free(cache_vars2.li);
     free(cache_vars2.fz);
     free(cache_vars2.zold);
     free(cache_vars2.DMin);
     free(cache_vars2.KK);
     free(cache_vars2.ux1);
     free(cache_vars2.ux2);

     return;

}




