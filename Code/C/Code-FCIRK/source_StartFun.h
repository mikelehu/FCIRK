/*-----------------------------------------------------------------------------*/
/*									       */
/*                        source_StartFun.h	         		       */
/*									       */
/* ----------------------------------------------------------------------------*/
{


/*------ Declarations --------------------------------------------------------*/

     Pkepler_sys *Pkepler;

/* ----- Implementation ------------------------------------------------------*/

     Pkepler=&params->Pkepler;

#if HIGH5 == 0
     KeplerFlowAll (neq,Pkepler->keplerkop,w,h/2,params);
#else
     KeplerFlowAll_high (neq,Pkepler->keplerkop,w,h/2,params);
#endif


#    ifdef DEBUGSTARTFUN

     printf("\nStartFun irteera neq=%i\n",neq);

     int i;
     double norma=0.,daux;

     norma=0.;
     for (i=0; i<neq-1;i++) norma+=w->uu[i]*w->uu[i];
     norma=sqrt(norma);

     printf("norma=%.20lg,",norma);
     printf("\n");


     printf("u ****\n");
     for (i=0; i<neq/2;i++)
     {
         daux=w->uu[i];
         printf("%.20lg,",daux);
     }

     printf("\n");

     for (i=neq/2; i<neq-1;i++)
     {
         daux=w->uu[i];
         printf("%.20lg,",daux);
     }

     printf("\n");


#    endif

     return;
}
