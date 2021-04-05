/*-----------------------------------------------------------------------------*/
/*									       */
/*                        source_OutputFun.h	         		       */
/*									       */
/* ----------------------------------------------------------------------------*/
{

/*------ Declarations --------------------------------------------------------*/

     int neq;
     neq=system->neq;

     Pkepler_sys *Pkepler;


/* ----- Implementation ------------------------------------------------------*/


     Pkepler=&params->Pkepler;

     if (((cache_stat->stepcount % options->sampling) == 0) || (cache_stat->laststep))
     {
          
            solution ux;
            ux.uu = (highprec *)malloc((neq)*sizeof(highprec)); 
            ux.uul = (val_type *)malloc((neq)*sizeof(val_type)); 
            double DH;
            int i;

            struct rec
            {
             __float128 t;
            __float128 uu[neq];
            };

           struct rec my_record;

           cache_stat->nout++;
  
           for(i=0; i<neq; i++)
           {
            ux.uu[i]=w->uu[i];
            ux.uul[i]=w->uul[i];
           }
      
           if(cache_stat->stepcount==0)
           {

               DH=0.;
               cache_stat->E0=system->ham(neq,&ux,params);
#ifdef DEBUG
               int n;
               int width = 46;
               char buf[128];
               printf("NS=%i,Initial energy: ",method->ns);
               n = quadmath_snprintf(buf, sizeof buf, "%+-#*.36Qe", width, cache_stat->E0);
               if ((size_t) n < sizeof buf) printf("%s\n",buf);
#endif
          }
          else
          {
#if HIGH == 0
               KeplerFlowAll (neq,Pkepler->keplerkop,&ux,-h/2,params);
#else
               KeplerFlowAll_high (neq,Pkepler->keplerkop,&ux,-h/2,params);
#endif
               DH=(system->ham(neq,&ux,params)-cache_stat->E0)/cache_stat->E0;
               if (FABS(DH)>cache_stat->MaxDE) cache_stat->MaxDE=FABS(DH);
          }


#ifdef DEBUG
         double aux;
         aux=t;
         printf("%i,%lg,%i, %lg\n",cache_stat->stepcount,aux,cache_stat->itcount,DH);
#endif

         my_record.t=t;
         for (i=0; i<neq;i++) my_record.uu[i]=ux.uu[i];
         fwrite(&my_record, sizeof(struct rec), 1, myfile);


#    ifdef DEBUGOUTPUTFUN

         double daux;
         printf("u=");
         for (i=0; i<neq/2; i++)
         {
           daux=ux.uul[i];
           printf("%.20lg,",daux);
         }
         printf("\n");
         for (i=neq/2; i<neq; i++) 
         {  
           daux=ux.uul[i];
           printf("%.20lg,",daux);
         }
         printf("\n");

#        endif


        free(ux.uu);
        free(ux.uul);

     }



     return;

}
