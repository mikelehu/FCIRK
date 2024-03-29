/*-----------------------------------------------------------------------------*/
/*									       */
/*                        source_FP_iteration.h	       		               */
/*									       */
/* ----------------------------------------------------------------------------*/

{

#ifdef PARALLEL
     int extern thread_count;
#endif

/* ----- First initializations -----------------------------------------------*/

     int neq,ns;
     BASE *z,*zold,*li,*fz,*DMin;
     BASE *ttau,nrm;

     neq=system->neq;
     ns=method->ns;

#if HIGH ==0 
     parameters *params;
#else
     parameters_high *params;
#endif

     z=cache_vars->z;
     zold=cache_vars->zold;
     li=cache_vars->li;
     fz=cache_vars->fz;
     DMin=cache_vars->DMin;
     ttau=method->ttau;
     params=&system->params;
     nrm=options->nrmdigits;

/*------ Declarations --------------------------------------------------------*/

     int i,is,js,isn,in,jsn,jsni,ism,ismi;
     val_type sum,dY;     
     int DD0,myfcn;
     bool cont,eval;
     bool plusIt;

/* ----- Implementation ------------------------------------------------------*/

 
     if (*D0<0) plusIt=false;
         else plusIt=true;

     DD0=0;         
     cont=false;      
     myfcn=0;


     for (is = 0; is<ns; is++)
     {
           isn=neq*is;

           for (i = 0; i<neq; i++)
           {
                 in=isn+i;
                 sum=u->ee[i]+method->m[ns*is]*li[i];
                 for (js =1; js<ns; js++)
                 {
                       jsn=ns*is+js;
                       jsni=neq*js+i;
                       sum+=method->m[jsn]*li[jsni];
                 }

                 zold[in]=z[in];
                 z[in]=u->uu[i]+sum;
           }
     }

     ism=ii*ns;
#ifdef PARALLEL
#      pragma omp parallel for num_threads(thread_count) \
                               default(shared)\
                               private(is,i,isn,in,dY,ismi) \
                               reduction (+:DD0,myfcn) \
                               reduction (||:cont)
#endif
     for (is = 0; is<ns; is++)
     {
           isn=neq*is;
           ismi=ism+is;
           eval=false;    

           for (i=0; i<neq; i++) 
           {   
               in=isn+i;
               
               /* Stop Criterion */
#if HIGH ==0                
	       dY=KFABS(Rmdigits(z[in],nrm)-Rmdigits(zold[in],nrm));
#else
	       dY=KFABS(Rmdigits_high(z[in],nrm)-Rmdigits_high(zold[in],nrm));
#endif	       

               if (dY>0.)
               {
                  eval=true;
                  if (dY<DMin[in])
                  {
                        DMin[in]=dY;
                        cont =true;
                  }
                  else if (dY>FABS(z[in])*1e-8)
                  {
                        cont =true;  
                  }
               }
               else
               {
                  DD0++;

//               2021-11-22 Rmdigits segurua izateko                              
	          if (KFABS(z[in]-zold[in])>0.) eval=true;	 
               
               }   

           }  
       
           if (eval)
           {
              myfcn++;
              system->f(neq,tn+method->hc[is],ttau[ismi],&z[isn],&fz[isn],params,cache_vars);

              for (i=0; i<neq; i++)
              {
                   in=isn+i;
                   li[in]=fz[in]*method->hb[is];
              }
           }
     }

     *D0=DD0;
     *iter0=cont;
     cache_stat->fcn+=myfcn;


     if (*iter0==false && *D0<(ns*neq) && plusIt)
     {
          *D0=-1;
          *iter0=true;
     }


#ifdef  NDEBUG


//  printf("FP_Iteration\n");
  printf("step=%i, iteration=%i, iter0=%i, D0=%i, plusIt=%i ****\n",
          cache_stat->stepcount, cache_stat->itcount, *iter0, DD0, plusIt); 
#endif

     
     return(0);


}


