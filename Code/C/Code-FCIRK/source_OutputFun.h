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
     
     __float128 Ei;
      


/* ----- Implementation ------------------------------------------------------*/


     Pkepler=&params->Pkepler;

     if (((cache_stat->stepcount % options->sampling) == 0) || (cache_stat->laststep))
     {
          
            solution ux;
            ux.uu=cache_vars->ux1;
            ux.ee=cache_vars->ux2;
            
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
            ux.ee[i]=w->ee[i];
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
               Ei=system->ham(neq,&ux,params);
               DH=(Ei-cache_stat->E0)/cache_stat->E0;
               if (FABS(DH)>cache_stat->MaxDE) cache_stat->MaxDE=FABS(DH);
          }


#ifdef DEBUG
         double aux;
         aux=t;
         printf("%i,%lg,%i, %lg\n",cache_stat->stepcount,aux,cache_stat->itcount,DH);
#endif

         my_record.t=t;
         for (i=0; i<neq;i++)
         {
          my_record.uu[i]=ux.uu[i];
          my_record.uu[i]+=ux.ee[i];
         
         }
         fwrite(&my_record, sizeof(struct rec), 1, myfile);       
 
        if (options->errorsL==true)
        {         
         my_record.t=t;
         for (i=0; i<neq;i++)
         {
          my_record.uu[i]=errj[i];
          errj[i]=0.;
         }
         fwrite(&my_record, sizeof(struct rec), 1, myfileER);
         }
         
      
         struct rec2
            {
              double t;
              double MaxRC;
              double MinRC;
              double Mean_RC;
              double Var_RC;
         };   
         
         
         struct rec2 my_record2;    
         
         if (cache_stat->stepcount>0)
         {
            my_record2.t=t;
            my_record2.MaxRC=cache_vars->MaxRC;
            my_record2.MinRC=cache_vars->MinRC;
            cache_vars->MaxRC=0.;
            cache_vars->MinRC=1.e6;
            my_record2.Mean_RC=cache_vars->RCSum/cache_stat->stepcount;
            my_record2.Var_RC=SQRT(cache_vars->RCSum2/cache_stat->stepcount-my_record2.Mean_RC*my_record2.Mean_RC);
        
            fwrite(&my_record2, sizeof(struct rec2), 1, myfileRC);
         }
 
   

#    ifdef DEBUGX

         double daux1,daux2,daux3,daux4;


         if (cache_stat->stepcount>0)
         {
          daux1=my_record2.t;
          daux2=my_record2.RC;
          daux3=my_record2.Mean_RC;
          daux4=my_record2.Var_RC;

          printf("RCFile.%i, %lg, %lg, %lg, %lg\n",cache_stat->stepcount,daux1,daux2,daux3,daux4);
         }
         
#        endif


     }



     return;

}
