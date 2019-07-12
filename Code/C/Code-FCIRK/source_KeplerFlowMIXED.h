/*-----------------------------------------------------------------------------*/
/*									       */
/*                        source_KeplerFlowMIXED.h         		       */
/*									       */
/* ----------------------------------------------------------------------------*/
{

/*------ Declarations --------------------------------------------------------*/

     int i;
     int dim=3;

     BASE3 r0,eta,alpha,beta,gamma,zeta;
     BASE3 r,rinv,X;
     BASE3 b[4],G[3];   
    
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

#if HIGH3 ==0

     r0=KSQRT(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]);
     eta=q[0]*v[0]+q[1]*v[1]+q[2]*v[2];
     alpha=k/r0;
     beta=2.*alpha-(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
     gamma=k/beta;
     zeta=k-beta*r0

     KeplerSolveE(r0,gamma,eta,beta,k,dt,&X,G);		
#else
     val_type kL,r0L,etaL,betaL,gammaL;
     val_type dtL,XL,GG[3];   

     r0=KSQRT(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]);  
     eta=q[0]*v[0]+q[1]*v[1]+q[2]*v[2];
     alpha=k/r0;
     beta=2.*alpha-(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
     gamma=k/beta;
     zeta=k-beta*r0;

     r0L=r0;
     gammaL=gamma;
     etaL=eta;
     betaL=beta;
     kL=k;
     dtL=dt;
     KeplerSolveE(r0L,gammaL,etaL,betaL,kL,dtL,&XL,GG);	

     X=XL;
     for (i=0; i<3; i++) G[i]=GG[i];
     KeplerSolveITER_high 
                 (r0,gamma,eta,beta,k,dt,&X,G); 

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

