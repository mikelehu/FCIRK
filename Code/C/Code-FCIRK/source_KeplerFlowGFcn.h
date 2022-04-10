/*-----------------------------------------------------------------------------*/
/*									       */
/*                        source_KeplerFlowGFcn.h   	       		       */
/*									       */
/* ----------------------------------------------------------------------------*/
{

     int dim=3;
     int nd=neq/2;

     int j,j1,j2;

     BASE u[neq],g[neq];


#if HIGH ==0 
     flowaux aux[keplerkop];
#else
     flowaux_high aux[keplerkop];
#endif


     for (j=0; j<keplerkop; j++)
     {
          j1=j*dim;
          j2=j1+nd;
#if HIGH ==0 
          KeplerFlowGen (dt, &U[j1], &U[j2], k[j], &u[j1], &u[j2], &aux[j]);
#else
          KeplerFlowGen_high (dt, &U[j1], &U[j2], k[j], &u[j1], &u[j2], &aux[j]);
#endif

     }

     GFcn (neq, t, u, g, params);
     for (j=0; j<neq; j++) G[j]=g[j];


     for (j=0; j<keplerkop; j++)
     {
          j1=j*dim;
          j2=j1+nd;
#if HIGH ==0 
          KeplerFlowGFcnAux (dt,&U[j1],&U[j2],&g[j1],&g[j2],&aux[j],k[j],&G[j1],&G[j2]);
#else
          KeplerFlowGFcnAux_high (dt,&U[j1],&U[j2],&g[j1],&g[j2],&aux[j],k[j],&G[j1],&G[j2]);
#endif

     }



#    ifdef DEBUGKEPLERFLOWGFCN

   printf("KeplerFlowGFcn, neq=%i\n",neq);
   
   int i;
   val_type norma;

   printf("dt=%lg\n",dt);

   norma=0.;
   for (i=0; i<neq-1; i++) norma+=G[i]*G[i];
   printf("Norm(f)=%lg\n",sqrt(norma));
   printf("\n");

#    endif

   return;

}





