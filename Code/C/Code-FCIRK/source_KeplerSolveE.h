/*-----------------------------------------------------------------------------*/
/*									       */
/*                        source_KeplerSolve.h	         		       */
/*									       */
/* ----------------------------------------------------------------------------*/
{
/*   gamma X + eta G2 + zeta G3 == dt -----> X */

/*

  Elliptic case (beta>0)
  F(X) := r0*G0 + eta0*G2 + k*G3 - t == 0 ->   X
  F'(X) =  r0*G0 + eta*G1 + k*G2
  |F''(X)| <= om := sqrt(eta^2 + (k -beta*r0)^2/beta) = k*e/sqrt(beta)
  An initial guess X0 \approx X is obtained as the (exact or approximate) zero
  of a cubic Taylor polynomial of F(X) centered at X = XX, where
    - XX = 0.  if  the Kantorovich’s alpha(0) <= 0.1*e
    - XX = 2*F(0)/F'(0)  if  alpha(0) <= 0.5*e
   - XX = (beta*t - eta)/k  otherwise

*/

/*------ Declarations --------------------------------------------------------*/


     int i;
     bool cont;

     BASE nu,om;
     BASE X0,x,dX,XX,dt,r;
     BASE G3,zeta0,eta,alpha,tau;
     BASE errX; 
     BASE S,S2,C,W;
     BASE a0,a1,a2,a3,a4,aux; 
     BASE b0,b1,b2, P0,Pd0;
     BASE c,c2,p,q,pb3,diskr;
     BASE lag1,lag2,lag3;
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
     zeta0 =k-beta*r0;
     om=KSQRT(eta0*eta0+zeta0*zeta0/beta);

#ifdef TRACE
     double e,E0,M0;
     e = om*nu/k;
     E0=KATAN2(nu*eta0,zeta0);
     M0 = E0 - e*KSIN(E0); 
     printf("M=%lg, E0=%lg, e=%lg \n", M0 + KPOW(nu,3)/k*t,E0,e);
#endif


     dX=t/r0;        // -F(0)/F'(0)
     tau = k/nu*FABS(dX/r0);


     if (tau<0.03)
     { 
         X0=0.; 
         a0=-t;
         a1=r0;
         a2=eta0;
         a3=zeta0;
         i=0;   
     } 
      else
     {
          if (tau<0.5) X0=dX;
          else X0=(beta*t-eta0)/k;

          x=nu*X0;
          G[0]=KCOS(x);
          G[1]=KSIN(x)/nu;
          G[2]=(1.-G[0])/beta;
          G3=(X0-G[1])/beta;
          a0=r0*G[1]+eta0*G[2]+k*G3-t;
          aux=r0*G[0]+eta0*G[1];
          a1=aux+k*G[2];
          a2=eta0*G[0]+zeta0*G[1];
          a3=-beta*aux+k*G[0];
          i=1;
     }




#ifdef TRACE
          printf("Taylor polynomial centered at EE0=%.20lg, X=%.20lg, tau=%.20lg, <=0.5 or <=0.03\n",
                  E0+nu*X0,X0,tau); 
#endif

     if ((KFABS(a1*beta)>=0.4*k) || tau<=0.5)
     { 
           a2*=0.5;
           a3*=0.16666666666666667;
           a4=-beta*a2/12.;
           dX=-a0/a1;     
           dX=-a0/(a1+a2*dX);    // One Halley's iteration
           b2=a3+dX*a4;
           b1=a2+dX*b2;
           b0=a1+dX*b1;
           P0=a0+dX*b0;
           Pd0=b0+dX*(b1+dX*(b2+dX*a4));
           dX-=P0/Pd0;           // Newton iteration
           X0+=dX;

#ifdef TRACE
           printf("Starting guess (approx.: X0 =%lg, dX=%lg, p(dX)=%lg\n", 
                    X0,dX,a0+a1*dX+a2*dX*dX+a3*dX*dX*dX+a4*dX*dX*dX*dX);
#endif
     }
      else
     { 
           if (a0==0 && a1==0)
           {  
#ifdef TRACE
               printf("Success: E=%lg, X=%lg, errX=%lg, iters=%i\n",
                       (E0 + KSQRT(beta)*X0),X0,0.,1);
#endif
               cont=false;
           }
           else
           {
               c=-a2/a3;
               c2=c*c;
               p=2.*a1/a3-c2;
               q=-3.*a0/a3-c*(1.5*p+0.5*c2); 
//             x^3 + 3*p*x - 2*q = 0  
               pb3=p*p*p;
               diskr=q*q+pb3;
               lag1=KSQRT(diskr);
               if (q>0)
                   lag2=q+lag1;
               else lag2=q-lag1;
               lag3=-pb3/lag2;
#if HIGH ==0  
               x=KPOW(lag2,1./3)+sign(lag3)*pow(KFABS(lag3),1./3);
#else
               x=KPOW(lag2,1./3)+sign_high(lag3)*pow(KFABS(lag3),1./3);
#endif
               X0+=c+x;
               X0+=dX;
#ifdef TRACE
               printf("Starting guess (Exact: X0 =%lg, dX=%lg, p(dX)=%lg\n ",
                       X0,dX,a0+a1*dX+a2/2*dX*dX+a3/(6*dX*dX*dX));
#endif         
          }
      }  


     
      x=0.5*nu*X0;
      KSINCOS(x,&S,&C);    
      S2=2.*S;
      W=S2*S;
      G[0]=1.-W;
      G[1]=S2*C/nu;
      G[2]=W/beta;
#ifdef TRACE
      G3=(X0-G[1])/beta;
      dt=t-(r0*G[1]+eta0*G[2]+k*G3);
      r=r0*G[0]+eta0*G[1]+k*G[2];
      dX=dt/r;
      alpha=om*FABS(dX/r);
      printf("alpha=%lg\n",alpha);
      printf("i=%i, E=%.16lg, X=%.16lg\n",
              i,E0+sqrt(beta)*X0,X0); 
#endif
    
//   At this point, we have good starting values for (X)  
     while (cont && i<imax)
     {
         G3=(X0-G[1])/beta;
         dt=t-(r0*G[1]+eta0*G[2]+k*G3);
         r=r0*G[0]+eta0*G[1]+k*G[2];
         dX=dt/r;
         i+=1;
         alpha=om*FABS(dX/r);
         errX=KFABS(0.52*alpha*dX);    // Error bound for Newton approximation,
                                       // under the assumption that alpha <= 0.032
         XX=X0+dX;
         Xmin=XX-errX;
         Xmax=XX+errX;

         if (Xmin==Xmax)
         {
           X0=XX;
           cont = false;
         } 
         else
         {
           eta=eta0*G[0]+zeta0*G[1]; 
           dX/=1.+0.5*dX*eta/r;
           X0+=dX; 
         }
#ifdef TRACE
         printf("i=%i, E=%lg, X=%lg, dE=%lg, alpha=%lg, errE=%lg\n",
                 i,E0+KSQRT(beta)*X0,X0, dX*nu, alpha, errX*nu); 
#endif

         x=0.5*nu*X0;     
         KSINCOS(x,&S,&C);    
         S2=2.*S;
         W=S2*S;
         G[0]=1.-W;
         G[1]=S2*C/nu;
         G[2]=W/beta; 

    }  


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

