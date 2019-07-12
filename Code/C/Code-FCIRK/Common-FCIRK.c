/*----------------------------------------------------------------------------*/
/*									      */
/*    Common-FCIRK.c							      */
/*									      */
/*	Functions: 							      */
/*	 InitStat()							      */
/*	 NormalizedDistance()						      */
/*	 Fixed_point_Step()						      */
/*	 General_FP_It()						      */
/*	 Summation()							      */
/*	 Main_FCIRK()							      */
/*       								      */
/* ---------------------------------------------------------------------------*/


#include <Common-FCIRK.h>
#include <float.h>


/******************************************************************************/
/* 					   				      */
/* InitStat 		         	  				      */
/*                                        				      */
/*									      */
/******************************************************************************/

void InitStat (ode_sys *system, gauss_method *gsmethod,
               solver_stat *thestatptr)
{

     int i,is,in,neq,ns;

     ns=gsmethod->ns;
     neq=system->neq;

     thestatptr->laststep = false;
     thestatptr->stepcount = 0;
     thestatptr->itcount=0;
     thestatptr->totitcount=0;
     thestatptr->totitcountzero=0;
     thestatptr->maxitcount=0;
     thestatptr->itzero=0;
     thestatptr->fcn=0;
     thestatptr->convergence=SUCCESS;
     thestatptr->nout=0;
     thestatptr->MaxDE=0.;

     thestatptr->z = (val_type *)malloc(neq*(ns)*sizeof(val_type));
     thestatptr->li = (val_type *)malloc(neq*ns*sizeof(val_type));
     thestatptr->fz = (val_type *)malloc(neq*ns*sizeof(val_type));
     thestatptr->zold = (val_type *)malloc(neq*ns*sizeof(val_type));

     for (is=0; is<ns; is++)
     {
          in=is*neq;
          for (i=0; i<neq; i++)
          {
               thestatptr->li[in+i]=0.;
               thestatptr->fz[in+i]=0.;
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


/*---------------- Declarations ----------------------------------------------*/

     int i,is,ix;
     val_type maxi,mayi,relerrors;
     val_type maxz,maxzold;

/* --------------- Implementation --------------------------------------------*/

     maxi=0.;

     for (i=0; i<neq; i++)
     {
          maxz=0.;
          maxzold=0.;

          for (is=0; is<ns;is++)
          {
               ix=neq*is+i;
               maxz=FMAX(maxz,FABS(z[ix]));
               maxzold=FMAX(maxzold,FABS(zold[ix]));
          }

          mayi=(maxz+maxzold)/2;
          relerrors=0.;

          for (is=0; is<ns; is++)
          {
               ix=neq*is+i;
               relerrors=FMAX(relerrors,
                              FABS(z[ix]-zold[ix])/
                                   (mayi*options->rtol[i]+options->atol[i]));
          }

          maxi=FMAX(maxi,relerrors);

     }

     return maxi;

}



/******************************************************************************/
/* 									      */
/*   Fixed_point_Step							      */
/* 									      */
/******************************************************************************/


void Fixed_point_Step      (  ode_sys *system,  solution *u,
                              val_type tn,  val_type h,
                              toptions *options,  gauss_method *method,
                              solver_stat *thestatptr)

{

#ifdef PARALLEL
     int extern thread_count;
#endif

/* ----- First initializations -----------------------------------------------*/
     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------ Declarations --------------------------------------------------------*/

     int i,is,isn,in;
     bool iter0;
     int D0;
     double difftest;
     val_type DMin[neq*ns];
     val_type *z,*zold,*fz,*li;
     val_type *ttau;     
     parameters *params;

/* ----- Implementation ------------------------------------------------------*/


     z=thestatptr->z;
     zold=thestatptr->zold;
     fz=thestatptr->fz;
     li=thestatptr->li;
     ttau=method->ttau;

     params=&system->params;

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

     thestatptr->itcount=0;

#ifdef PARALLEL
#      pragma omp parallel for num_threads(thread_count) private(isn)
#endif
     for (is = 0; is<ns; is++)
     {
           isn=neq*is;
           thestatptr->fcn++;         
           system->f(neq,tn+method->hc[is],ttau[is],&z[isn],&fz[isn],params);
           for (i=0; i<neq; i++) li[isn+i]=fz[isn+i]*method->hb[is];
     }

     thestatptr->itcount++;

     while (iter0 && thestatptr->itcount<MAXIT)
     {
           General_FP_It (system,u,tn,h,method,thestatptr,&D0, &iter0, DMin);      
           thestatptr->itcount++;

     }


     if (thestatptr->itcount==MAXIT)
     {
           printf("Break: step(MAXIT)=%i\n",thestatptr->itcount);
           thestatptr->convergence=FAIL;
     }

     if (D0<(ns*neq))
     {
           difftest=NormalizedDistance(neq,ns,options,z,zold);
           if (difftest>1.)
           {
                thestatptr->convergence=FAIL;
                printf("Lack of convegence of Fixed point iteration:\
                        step=%i,iteration=%i,",
                        thestatptr->stepcount,thestatptr->itcount);
                printf("difftest=%lg\n",difftest);
           }
     }
     else
     {
           (thestatptr->totitcountzero)+=(thestatptr->itcount);
           thestatptr->itzero++;
     }

     return;

}



/******************************************************************************/
/*									      */
/*      General_FP_It: General Fixed Point Iteration			      */
/*									      */
/******************************************************************************/


int General_FP_It      ( ode_sys *system,  solution *u,  val_type tn,
                         val_type h,  gauss_method *method,
                         solver_stat *thestatptr,
                         int *D0, bool *iter0, val_type *DMin)
{

#ifdef PARALLEL
     int extern thread_count;
#endif

/* ----- First initializations -----------------------------------------------*/

     int neq,ns;
     parameters *params;
     val_type *z,*zold,*li,*fz;
     val_type *ttau;

     neq=system->neq;
     ns=method->ns;
     params=&system->params;

     z=thestatptr->z;
     zold=thestatptr->zold;
     li=thestatptr->li;
     fz=thestatptr->fz;
     ttau=method->ttau;

/*------ Declarations --------------------------------------------------------*/

     int i,is,js,isn,in,jsn,jsni;
     val_type sum,dY;     
     int DD0,myfcn;
     bool cont,eval;
     bool plusIt;

/* ----- Implementation ------------------------------------------------------*/

 
     if (*D0<0) plusIt=false;
         else plusIt=true;

     DD0=0;         
     cont=false; 
     eval=false;       
     myfcn=0;


     for (is = 0; is<ns; is++)
     {
           isn=neq*is;

           for (i = 0; i<neq; i++)
           {
                 in=isn+i;
                 sum=method->m[ns*is]*li[i];
                 for (js =1; js<ns; js++)
                 {
                       jsn=ns*is+js;
                       jsni=neq*js+i;
                       sum+=method->m[jsn]*li[jsni];
                 }

                 zold[in]=z[in];
                 z[in]=u->uul[i]+sum;
           }
     }


#ifdef PARALLEL
#      pragma omp parallel for num_threads(thread_count) \
                               default(shared)\
                               private(is,i,isn,in,dY) \
                               reduction (+:DD0,myfcn) \
                               reduction (||:cont)
#endif
     for (is = 0; is<ns; is++)
     {
           isn=neq*is;
           for (i=0; i<neq; i++) 
           {   
               in=isn+i;
               
               /* Stop Criterion */
	       dY=FABS(z[in]-zold[in]);

               if (dY>0.)
               {
                  eval=true;
                  if (dY<DMin[in])
                  {
                        DMin[in]=dY;
                        cont =true;
                  }
               }
               else
               {
                  DD0++;
               }   

           }  
       
           if (eval)
           {
              myfcn++;
              system->f(neq,tn+method->hc[is],ttau[is],&z[isn],&fz[isn],params);

              for (i=0; i<neq; i++)
              {
                   in=isn+i;
                   li[in]=fz[in]*method->hb[is];
              }
           }
     }

     *D0=DD0;
     *iter0=cont;
     thestatptr->fcn+=myfcn;


     if (*iter0==false && *D0<(ns*neq) && plusIt)
     {
          *D0=-1;
          *iter0=true;
     }


#ifdef MDEBUG


  printf("General_FP_IIT\n");
  printf("step=%i, iteration=%i, iter0=%i, D0=%i, plusIt=%i ****\n",
          thestatptr->stepcount, thestatptr->itcount, *iter0, DD0, plusIt); 
#endif

     
     return(0);


}



/******************************************************************************/
/* 									      */
/*    Summation								      */
/* 									      */
/* 									      */
/******************************************************************************/

void Summation                ( gauss_method *gsmethod,
                                solution *u, ode_sys *system,
                                toptions *options, solver_stat *thestatptr)

{


/* ----- First initializations -----------------------------------------------*/

     int neq,ns;

     neq=system->neq;
     ns=gsmethod->ns;

/*------ Declarations --------------------------------------------------------*/

     int i,is,isn;
     val_type *fz;
     fz=thestatptr->fz;


/* ----- Implementation ------------------------------------------------------*/


//---- High-prec computation of Li 


    for (is=0; is<ns; is++)
    {
        isn=neq*is;
        u->uu[0]+=fz[isn]*gsmethod->hb[is];

        for (i = 1; i<neq; i++)
                u->uu[i]+=fz[isn+i]*gsmethod->hb[is];
    }


    for (i = 0; i<neq; i++) u->uul[i]=u->uu[i];


#ifdef MDEBUG

   double aux;

   aux=0.;
   for (i=0; i<neq; i++) aux+=u->uu[i]*u->uu[i];
   printf("Norm(w)=%lg\n",sqrt(aux));

#endif


    return;

}




/******************************************************************************/
/* 									      */
/*   Main-FCIRK	  							      */
/*   									      */
/* 									      */
/******************************************************************************/

void Main_FCIRK
(val_type t0, val_type t1, val_type h,
  gauss_method *gsmethod, solution *u,
  ode_sys *system, toptions *options,
  solver_stat *thestatptr)

{

/* ----- First initializations -----------------------------------------------*/

     int neq;

     Pkepler_sys *Pkepler;
     parameters *params;

     neq=system->neq;
     params=&system->params;
     Pkepler=&params->Pkepler;

/*------ Declarations --------------------------------------------------------*/

     FILE *myfile;

     int i,is,in,ns,istep,nstep;
     int js,jsn;
     val_type tn;
     val_type *z,*coef;
     val_type *zold,*fz;
     val_type sum;

     fz=thestatptr->fz;
     z=thestatptr->z;
     zold=thestatptr->zold;
     ns=gsmethod->ns;
     coef=gsmethod->nu;

/* ----- Implementation ------------------------------------------------------*/

#ifdef IOUT
     myfile = fopen(options->filename,"wb");
#endif

     tn=t0;
     nstep=round((t1-t0)/h);
     thestatptr->execution=SUCCESS;

#ifdef IOUT
     options->TheOutput(system,gsmethod,tn,h,u,thestatptr,params,options,myfile);
#endif

     system->StartFun(neq,t0,h,u,params);

     for(istep=0; (istep<nstep); istep++)
     {
 

          if (istep==0 || h>7)

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


          Fixed_point_Step (system,u,tn,h,options,gsmethod,thestatptr);
          if (thestatptr->convergence==FAIL || isnan(fz[0]))
          {
               thestatptr->execution=FAIL;
               nstep=istep;

          }
          else
          {
               Summation (gsmethod,u,system,options,thestatptr);
               tn=t0+(istep+1)*h;

               thestatptr->stepcount++;
               (thestatptr->totitcount)+=(thestatptr->itcount);
               if ((thestatptr->itcount)>(thestatptr->maxitcount))
                  (thestatptr->maxitcount)=(thestatptr->itcount);

               system->ProjFun(neq, tn, h, u, params);
#ifdef IOUT
               options->TheOutput(system,gsmethod,tn,h,u,thestatptr,params,options,myfile);
#endif
               if (tn+h>=t1 && thestatptr->laststep==false)
                    thestatptr->laststep=true;

          }
     }

#ifdef DEBUG
     if (thestatptr->execution==SUCCESS) printf("Success execution!!\n");
     else printf("Execution Fail!!\n");
#endif

#ifdef IOUT
     fclose(myfile);
#endif

     return;

}




