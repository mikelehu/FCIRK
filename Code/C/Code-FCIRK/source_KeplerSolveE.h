/*-----------------------------------------------------------------------------*/
/*									          */
/*                        source_KeplerSolve.h	                         */
/*   2021-06-18                                                                */
/*									          */
/* ----------------------------------------------------------------------------*/
{

/*------ Declarations --------------------------------------------------------*/


     bool cont;

     BASE nu,nuinv,om,bom,betainv;
     BASE X0,x,dX,discr;
     BASE G3,alpha0,delta0,deltaX,lambda0,lambdaX;
     BASE errX; 
     BASE C,S,S2;
     BASE a0,a1,a2,a3; 
     BASE c,c2,p,q,pb3,diskr;
     BASE aux,lag1,lag2;
     BASE F0,F0d,F0dd;
     BASE Xmin,Xmax;

     int imax=100;

#ifdef MDEBUG

    double rx0,gammax0,etax0,betax0,kx0,tx0;

    rx0=r0;
    gammax0=gamma0;
    etax0=eta0;
    betax0=beta;
    kx0=k;
    tx0=t;

#endif

/* ----- Implementation ------------------------------------------------------*/

#ifdef MDEBUG
    double mr,mgamma,meta,mbeta,mk,mt;
    mr=r0;
    mgamma=gamma0;
    meta=eta0;
    mbeta=beta;
    mk=k;
    mt=t;
    printf("keplerSolve input ******\n");
    printf("rx0=%lg,gammax0=%lg,etax0=%lg,betax0=%lg, kx0=%lg, tx0=%lg\n",mr,mgamma,meta,mbeta,mk,mt);

#endif
   
     cont=true;


     nu=KSQRT(beta);
     nuinv=1/nu;
     betainv=1/beta;
     bom=KSQRT(beta*eta*eta+zeta*zeta);  // sqrt(beta)*Omega
     om= bom*nuinv;   // Omega
     *iter=0;
     delta0 = -t/r0;
     alpha0 =om*KFABS(delta0/r0);
     
//    Computation of initial guess X

     if (alpha0<=0.3 && nu*KFABS(delta0) <= 0.4)
     { 
         aux = eta/r0;
         lambda0 = aux*delta0;
 //       One super-Halley iteration with X0=0
         X0= -(1 + lambda0/(2-2*lambda0))*delta0;
     } 
      else // Computation of alternative initial value X
     {
          X0 = (beta*t-eta)/k;
          x = nu*X0;         
          if (KFABS(x)>0.78)
          {
            KSINCOS(x,&S,&G[0]);         
            G[1] = S*nuinv;
            G[2]= (1 - G[0])*betainv;
          }
          else
          {
            KSINCOS(x/2,&S,&C);
            S2 = S*S;
            G[0] = 1 - 2*S2;
            G[1] = 2*S*C*nuinv;
            G[2]= 2*S2*betainv;
          }

          G3 = (X0-G[1])*betainv;
          a0 = r0*X0 + eta*G[2] + zeta*G3 - t;
          a1 = r0 + eta*G[1] + zeta*G[2];
          a2 = eta*G[0] + zeta*G[1];
          
          if (KFABS(a1 * beta)>=0.63*k)
          { 
/*
            Then, the solution of 
            if abs(a1 * beta) >= 0.62546*k #Then, the solution of 
            the quadratic equation a0 + a1* dX + a2/2 * dX^2 
            that is closest to dX = -a0/a1  is computed
*/        
            discr = a1*a1-2*a0*a2;
            dX = -2*a0/(a1 + KSQRT(discr));  
#ifdef TRACE
            printf("p2(dX)=%.16lg\n", a0 + a1* dX + a2/2 * dX*dX);
#endif        
            X0+=dX; 
          }
          else
          {
/*          Otherwise, the cubic equation 
            a0 + a1* dX + a2/2 * dX^2 + a3/6 * dX^3 = 0 is solved
*/            
            a3 = -beta*eta*G[1] + zeta*G[0];
            c = -a2/a3;
            c2 = c*c;
            p = 2*a1/a3 - c2;
            q = -3*a0/a3 - c/2 * (3*p + c2);
//           z^3 + 3*p*z - 2*q = 0;
            pb3 = p*p*p;
            diskr = q*q + pb3;
            lag1 = KSQRT(diskr);            
            if (q>0)   //cbrt
                lag2=KPOW(q+lag1,1./3);
            else lag2=KPOW(q-lag1,1./3);   
            dX=c+lag2-p/lag2;	  
 #ifdef TRACE
            printf("p3(dX)=%.16lg\n", a0 + a1 * dX + a2/2 * dX*dX + a3/6 * dX*dX*dX );
#endif                       
            X0+=dX;                
          }       
        *iter+=1; 
     }
//       End of computation of initial guess
                         
     x=nu*X0;
     if (KFABS(x)>=0.78)
     { 
         KSINCOS(x,&S,&G[0]);
         G[1] = S*nuinv;
         G[2] = (1-G[0])*betainv;
     }
     else
     {
         KSINCOS(x/2,&S,&C);
         S2 = S*S;
         G[0] = 1 - 2*S2;
         G[1] = 2*S*C*nuinv;
         G[2] = 2*S2*betainv;
     }   
             
     G3 = (X0-G[1])*betainv;
     cont = true;
     
     while (cont)
     {
         *iter+=1;
         F0 = r0*X0 + eta*G[2] + zeta*G3 - t;
         F0d = r0 + eta*G[1] + zeta*G[2];
         F0dd = eta*G[0] + zeta*G[1];
         deltaX = F0/F0d;
         aux = F0dd/F0d;
         lambdaX = aux*deltaX;
         dX = -(1 + lambdaX/(2-2*lambdaX))*deltaX;
         X0+= dX;   
         errX = bom*KFABS(dX*dX*dX/F0d) + 3*aux*(dX+deltaX)*(dX+deltaX);
         Xmin=X0-errX;
         Xmax=X0+errX;
#ifdef TRACE
         printf("iter=%i, X=%.16lg, dX=%.16lg, errX=%.16lg, Xmin=%.16lg, Xmax=%.16lg\n", *iter,X0,dX,errX,Xmin,Xmax);
#endif
         x=nu*X0;
         if (KFABS(x)>0.78)
         {
            KSINCOS(x,&S,&G[0]);          
            G[1] = S*nuinv;
            G[2] = (1 - G[0])*betainv;
         }
         else
         {
            KSINCOS(x/2,&S,&C);
            S2 = S*S;
            G[0] = 1 - 2*S2;
            G[1] = 2*S*C*nuinv;
            G[2] = 2*S2*betainv;
         }
         
         G3 = (X0 - G[1])*betainv;
         cont = (Xmin!=Xmax) && (*iter < imax);

     }  

 //    r = r0 + eta*G[1] + zeta*G[2];   ez da erabiltzen
     *X=X0;
    
    
#ifdef MDEBUG
    double rx,gammax,etax,betax,kx,tx;
    if (i==imax)
    {
        rx=r0;
        gammax=gamma0;
        etax=eta0;
        betax=beta;
        kx=k;
        tx=t;
        printf("keplerSolve convergence fail ******\n");
        printf("r0=%lg,gamma0=%lg,eta0=%lg,beta=%lg, k=%lg, t=%lg\n",rx,gammax,etax,betax,kx,tx);
    }
    else 
    {
       double G0,G1,G2,eX;
       G0=G[0];
       G1=G[1];
       G2=G[2];
       eX=*X;
       tx=t;
       printf("keplerSolve  convergence success, t=%lg\n",tx);
       printf("G0=%lg, G1=%lg, G2=%lg, X=%lg\n",G0,G1,G2,eX);
    }
#endif

    return;


}

