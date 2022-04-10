/*-----------------------------------------------------------------------------*/
/*									       */
/*                        source_KeplerFlow.h	         		       */
/*									       */
/* ----------------------------------------------------------------------------*/
{


/*------ Declarations --------------------------------------------------------*/

     int i, iter;
     int dim=3;

     BASE r0,eta,v02,alpha,beta;
     BASE zeta;           
     BASE r,rinv,X;
     BASE b[4],G[3];   
    
/* ----- Implementation ------------------------------------------------------*/
   

     if (k==0)				
     {

       for (i=0; i<dim; i++)
       {
          dq[i]=dt*v[i];
          dv[i]=0.;
       }

     }
     else
     {
       r0=KSQRT(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]);
       eta=q[0]*v[0]+q[1]*v[1]+q[2]*v[2];
       v02=(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
       alpha=k/r0;
       beta=2.*alpha-v02;
//       gamma=k/beta;
       zeta=k-beta*r0;

  
#if HIGH ==0
       KeplerSolveE(r0,eta,zeta,beta,k,dt,&X,G,&iter);     
#else
       KeplerSolveE_high(r0,eta,zeta,beta,k,dt,&X,G,&iter); 
#endif
       r= r0+eta*G[1]+zeta*G[2];   
       rinv=1./r;
     
       b[0]=-alpha*G[2];
       b[1]=r0*G[1]+eta*G[2];
       b[2]=-alpha*G[1]*rinv;
       b[3]=-k*G[2]*rinv;

       for (i=0; i<dim; i++)
       {
         dq[i]=b[0]*q[i]+b[1]*v[i];
         dv[i]=b[3]*v[i]+b[2]*q[i];
       }

    }

}

