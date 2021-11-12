/*-----------------------------------------------------------------------------*/
/*									       */
/*                        source_IRKstep_fixed.h	       	       */
/*									       */
/* ----------------------------------------------------------------------------*/

{
#ifdef PARALLEL
     int extern thread_count;
#endif

/* ----- First initializations -----------------------------------------------*/
     int neq,ns;

     neq=system->neq;
     ns=method->ns;

/*------ Declarations --------------------------------------------------------*/

     int i,is,isn,in,ism,ismi;
     bool iter0;
     int D0;
     double difftest;
     BASE *z,*zold,*fz,*li,*DMin;
     BASE *ttau; 


#if HIGH ==0 
     int jsn,js;
     parameters *params;
     Pkepler_sys *Pkepler;
     BASE sum; 
     BASE *coef;
     coef=method->nu;  
     params=&system->params;
     Pkepler=&params->Pkepler;
#else
     parameters_high *params;
     params=&system->params;
#endif

/* ----- Implementation ------------------------------------------------------*/




     z=cache_vars->z;
     zold=cache_vars->zold;
     fz=cache_vars->fz;
     li=cache_vars->li;
     DMin=cache_vars->DMin;
     ttau=method->ttau;    


/* --- Interpolation ---------------------------------------------------------*/

    if (cache_stat->interpolate !=true)
    {

        for (is = 0; is<ns; is++)
        {
            in=is*neq;
            for (i = 0; i<neq; i++) z[in+i]=u->uul[i];
        }
    }

    else
    {
#if HIGH ==0 
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
#else
      printf("Error. Interpolate=true in high precision\n");
      for (is = 0; is<ns; is++)
      {
            in=is*neq;
            for (i = 0; i<neq; i++) z[in+i]=u->uul[i];
      }
#endif
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
           system->f(neq,tn+method->hc[is],ttau[ismi],&z[isn],&fz[isn],params);
           for (i=0; i<neq; i++) li[isn+i]=fz[isn+i]*method->hb[is];
     }    

     cache_stat->itcount++;

     while (iter0 && cache_stat->itcount<MAXIT)
     {
#if HIGH ==0 
           FP_Iteration (system,u,tn,ii,h,options,method,cache_stat,cache_vars,&D0, &iter0); 
#else
           FP_Iteration_high (system,u,tn,ii,h,options,method,cache_stat,cache_vars,&D0, &iter0);
#endif     
           cache_stat->itcount++;

     }


     if (cache_stat->itcount==MAXIT)
     {
           printf("Break: step(MAXIT)=%i\n",cache_stat->itcount);
           cache_stat->convergence=FAIL;
     }

     if (D0<(ns*neq))
     {
#if HIGH ==0
           difftest=NormalizedDistance(neq,ns,options,z,zold);
#else
           difftest=NormalizedDistance_high(neq,ns,options,z,zold);
#endif 
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
#if HIGH ==0
         Summation (method,u,system,options,cache_vars);
#else
         Summation_high (method,u,system,options,cache_vars);
#endif 

     return;

}
