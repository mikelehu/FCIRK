/*-----------------------------------------------------------------------------*/
/*									       */
/*                        source_NbodyOde.h      	       		       */
/*									       */
/* ----------------------------------------------------------------------------*/
{


/*------ Declarations --------------------------------------------------------*/

   BASE *k;
#if HIGH ==0 
   Pkepler_sys *Pkepler;
#else
   Pkepler_sys_high *Pkepler;
#endif

/* ----- Implementation ------------------------------------------------------*/

   Pkepler=&params->Pkepler;
   k=Pkepler->K;   


#if HIGH ==0 
   KeplerFlowGFcn (NbodyGFcn,neq,t,u,Pkepler->keplerkop,k,params,ttau,f);   
#else
   KeplerFlowGFcn_high (NbodyGFcn_high,neq,t,u,Pkepler->keplerkop,k,params,ttau,f);  
#endif


#ifdef NDEBUG

   printf("\nNbodyODe, neq=%i,keplerkop=%i\n",neq,Pkepler->keplerkop);

   int i;
   double aux;
   aux=ttau;
   printf("ttau=%lg\n",aux);

   aux=0.;
   for (i=0; i<neq-1; i++) aux+=u[i]*u[i];
   printf("Norm(u)=%lg\n",sqrt(aux));
   aux=0;
   for (i=0; i<neq-1; i++) aux+=f[i]*f[i];
   printf("Norm(f)=%lg\n",sqrt(aux));
   printf("\n");

#  endif

   return;

};


