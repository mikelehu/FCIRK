/*-----------------------------------------------------------------------------*/
/*									       */
/*                        source_KeplerFlowGen.h   	       		       */
/*									       */
/* ----------------------------------------------------------------------------*/
{

     int id, iter;
     int dim=3;   


#ifdef XDEBUG         
     printf("KeplerFlowGen.....");    
               
     double q0,q1,q2;
     double v0,v1,v2;
     q0=Q[0]; q1=Q[1]; q2=Q[2];
     v0=V[0]; v1=V[1]; v2=V[2];
     printf("KeplerFlowGen   Q=(%lg,%lg,%lg), V=(%lg,%lg,%lg) !!!\n",q0,q1,q2,v0,v1,v2);

#endif


     if (k==0)
     {

        for (id=0; id<dim; id++)
        {
          qq[id]=Q[id]+dt*V[id];
          vv[id]=V[id];
        }

     }
     else
     {    
          
        aux->r0=KSQRT(Q[0]*Q[0]+Q[1]*Q[1]+Q[2]*Q[2]);
        aux->eta=Q[0]*V[0]+Q[1]*V[1]+Q[2]*V[2];
        aux->alpha=k/(aux->r0);
        aux->beta=2*(aux->alpha)-(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
        aux->gamma=k/(aux->beta);
        aux->zeta=k-(aux->beta)*(aux->r0);

#if HIGH ==0 
        KeplerSolveE(aux->r0,aux->eta,aux->zeta,aux->beta,k,dt,&aux->X,aux->GG,&iter);	
#else
        KeplerSolveE_high(aux->r0,aux->eta,aux->zeta,aux->beta,k,dt,&aux->X,aux->GG,&iter);
#endif
        aux->r= (aux->r0)+(aux->eta)*aux->GG[1]+(aux->zeta)*aux->GG[2];	        
        aux->rinv=1/(aux->r);
     
        aux->b[0]=-(aux->alpha)*aux->GG[2];
        aux->b[1]=(aux->r0)*aux->GG[1]+(aux->eta)*aux->GG[2];
        aux->b[2]=-(aux->alpha)*aux->GG[1]*(aux->rinv);
        aux->b[3]=-k*aux->GG[2]*(aux->rinv);

        for (id=0; id<dim; id++)
        {
           qq[id]=(1+aux->b[0])*Q[id]+aux->b[1]*V[id];
           vv[id]=(1+aux->b[3])*V[id]+aux->b[2]*Q[id];
        }

     }


#    ifdef NDEBUG

     int i;
     double xx;

     printf("KeplerFlowGen*******************\n");


//     printf("dt=%lg\n",dt);

//     printf("r0=%lg,eta=%lg,alpha=%lg,beta=%lg,gamma=%lg,zeta=%lg\n",
//             aux->r0,aux->eta,aux->alpha,aux->beta,aux->gamma,aux->zeta);
//     printf("X=%lg,G0=%lg, G1=%lg, G2=%lg, r=%lg, rinv=%lg\n",
//             aux->X,aux->GG[0],aux->GG[1],aux->GG[2],aux->r,aux->rinv);
//     printf("b11=%lg, b12=%lg, b21=%lg, b22=%lg\n", 
//             aux->b[0],aux->b[1],aux->b[2],aux->b[3]);

     printf("q,v=");
     for (i=0; i<dim; i++)
     {  xx=qq[i]; 
        printf("%lg,",xx);
     }
     printf("\n");
     for (i=0; i<dim; i++)
     {
        xx=vv[i];
        printf("%lg,",xx);
     }
     printf("\n");

#    endif

     return;

}          


