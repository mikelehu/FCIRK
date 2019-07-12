/*-----------------------------------------------------------------------------*/
/*									       */
/*                        source_ProjFun.h	         		       */
/*									       */
/* ----------------------------------------------------------------------------*/
{

/*------ Declarations --------------------------------------------------------*/

     Pkepler_sys *Pkepler;

/* ----- Implementation ------------------------------------------------------*/


     Pkepler=&params->Pkepler;

#if HIGH1 == 0
     KeplerFlowAll (neq,Pkepler->keplerkop,w,h,params);
#else
     KeplerFlowAll_high (neq,Pkepler->keplerkop,w,h,params);
#endif



#    ifdef DEBUGPROJFUN

     int i;
     val_type norma;

     norma=0.;
     for (i=0; i<neq-1; i++) norma+=w->uu[i]*w->uu[i];
     printf("Norm(w)=%Lg\n",sqrt(norma));
     printf("\n");

     printf("u ****\n");
     for (i=0; i<neq-1;i++) printf("%.4Lg,",w->uu[i]);
     printf("\ne ****\n");
     for (i=0; i<neq-1;i++) printf("%.4Lg,",w->ee[i]);
     printf("\n");

#    endif 


     return;
}
