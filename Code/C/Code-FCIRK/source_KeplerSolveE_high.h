/*-----------------------------------------------------------------------------*/
/*									          */
/*                        source_KeplerSolve_high.h	                         */
/*   2021-06-18                                                                */
/*									          */
/* ----------------------------------------------------------------------------*/
{

/*------ Declarations --------------------------------------------------------*/

     LOW r0_, eta_, zeta_, beta_, k_, t_;
     LOW X_, G_[3];
  
     BASE nu,nuinv,betainv;
     BASE X0,x,dX;
     BASE G3,deltaX,lambdaX;
     BASE C,S,S2;
     BASE aux;
     BASE F0,F0d,F0dd;


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

    r0_=r0;
    eta_=eta;
    zeta_=zeta;
    beta_=beta;
    k_=k;
    t_=t;
        
    KeplerSolveE (r0_, eta_, zeta_, beta_, k_ ,t_, &X_, G_, iter);
    

    X0=X_;     
    nu=KSQRT(beta);
    nuinv=1/nu;
    betainv=1/beta;                              
    x=nu*X0;
     
    if (KFABS(x)>=0.78q)
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
     *iter+=1;
     F0 = r0*X0+eta*G[2]+zeta*G3-t;
     F0d = r0+eta*G[1]+zeta*G[2];
     F0dd = eta*G[0]+zeta*G[1];
     deltaX = F0/F0d;
     aux = F0dd/F0d;
     lambdaX = aux*deltaX;
     dX = -(1 + lambdaX/(2-2*lambdaX))*deltaX;
     X0+= dX;   
#ifdef TRACE
     double aux1, aux2;
     aux1=X0;
     aux2=dX;
     printf("iter=%i, X=%.16lg, dX=%.16lg\n", *iter,aux1,aux2);
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
  
    *X=X0;
  
    return;

}

