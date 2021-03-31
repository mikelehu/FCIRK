/*-----------------------------------------------------------------------------*/
/*									       */
/*                        Julia_Source_FCIRK.h	         		       */
/*									       */
/* ----------------------------------------------------------------------------*/
{

/*------ initial declarations--------------------------------------------------*/

    int neq,nbody;
    int dim=3;

    neq=Nkepler*2*dim+Moreq;
    nbody=neq/(2*dim);

#ifdef DEBUG    
    printf("adaptive=%B,neq=%i, Nkekepler=%i, Moreq=%i\n",adaptive,neq,Nkepler,Moreq);
#endif
 
/*------ declarations --------------------------------------------------------*/

    int i;

    tmethod method;
    tcoeffs gsmethod;
    tcoeffs_h gsmethod_high;
    solution u;
    toptions options;
    tode_sys ode_system;
    ode_sys system;
    ode_sys_high system_high;
    Pkepler_sys Pkepler;
    Pkepler_sys_high Pkepler_high;
    tcache cache;
    tcache_stat cache_stat;
    tcache_vars cache_vars;
    tcache_vars_high cache_vars_high;

    clock_t clock0, clock1;
    time_t  wtime0,wtime1;

    val_type hq,t0q,t1q;



/* ----------- implementation  --------------------------------------------*/

    mpfr_set_default_prec (256);  // BigFloat's default precision in Julia


/* ----------- Initial values--------------------------------------------*/


    system.neq=neq;
    system_high.neq=neq;
    gsmethod.ns=ns;
    gsmethod_high.ns=ns;

    u.uu = (highprec *)malloc(neq*sizeof(highprec));
    u.uul = (val_type *)malloc(neq*sizeof(val_type));

    for (i=0; i<neq; i++)
    {
       u.uu[i]=mpfr_get_float128(*(u0+i),MPFR_RNDN);      
       u.uul[i]=u.uu[i];
    }

    hq=h;
    t0q=t0;
    t1q=t1;


#ifdef INITVALUES
   printf("Initial values: position\n");
   int id;
   for (i=0; i<nbody; i++)
   {
       for (id=0; id<dim; id++)
       {
          n = quadmath_snprintf(buf, sizeof buf, "%+-#*.36Qe", width, u.uu[i*dim+id]);
          if ((size_t) n < sizeof buf) printf("%s,",buf);

       }
       printf("\n");

   }

   printf("Initial values: velocity\n");
   for (i=0; i<nbody; i++)
   {
       for (id=0; id<dim; id++)
       {
          n = quadmath_snprintf(buf, sizeof buf, "%+-#*.36Qe", width, u.uu[nbody*dim+i*dim+id]);
          if ((size_t) n < sizeof buf) printf("%s,",buf);

       }
       printf("\n");

   }
#endif


/* -----------Parameters of the system -----------------------------------*/

    system.params.rpar =(val_type *)malloc(rlen*sizeof(val_type));
    system.params.rparhigh =(__float128 *)malloc(rlen*sizeof(__float128));
    system_high.params.rpar =(__float128 *)malloc(rlen*sizeof(__float128));
    system.params.ipar =(int *)malloc(ilen*sizeof(int));
    system_high.params.ipar =(int *)malloc(ilen*sizeof(int));
    system.params.numrpar=rlen;
    system_high.params.numrpar=rlen;
    system.params.numipar=ilen;
    system_high.params.numipar=ilen;

    for (i=0; i<rlen; i++)
    {
         system.params.rparhigh[i]=mpfr_get_float128(*(rpar+i),MPFR_RNDN);
         system_high.params.rpar[i]=mpfr_get_float128(*(rpar+i),MPFR_RNDN);
         system.params.rpar[i]=system.params.rparhigh[i];
    }


#ifdef INITVALUES

     int n;
     int width = 40;  
     char buf[128];

     printf("Rpar\n");
     for(i=0; i<rlen; i++)
     {
       n = quadmath_snprintf(buf, sizeof buf, "%+-#*.36Qe", width, system.params.rparhigh[i]);
       if ((size_t) n < sizeof buf) printf("%s\n",buf);
     }

     printf("Rpar (system_high\n");
     for(i=0; i<rlen; i++)
     {
       n = quadmath_snprintf(buf, sizeof buf, "%+-#*.36Qe", width, system_high.params.rpar[i]);
       if ((size_t) n < sizeof buf) printf("%s\n",buf);
     }
#endif

 
    for (i=0; i<ilen; i++)
    { 
       system.params.ipar[i]=ipar[i]; 
       system_high.params.ipar[i]=ipar[i];  
    }       

    Pkepler.K =(val_type *)malloc(nbody*sizeof(val_type));
    Pkepler.K2 =(val_type *)malloc(nbody*sizeof(val_type));
    Pkepler.Khigh =(__float128 *)malloc(nbody*sizeof(__float128));
    Pkepler.Mu =(val_type *)malloc(nbody*sizeof(val_type));
    Pkepler.Muhigh =(__float128 *)malloc(nbody*sizeof(__float128));

    Pkepler_high.K =(__float128 *)malloc(nbody*sizeof(__float128));
    Pkepler_high.K2 =(__float128 *)malloc(nbody*sizeof(__float128));
    Pkepler_high.Mu =(__float128 *)malloc(nbody*sizeof(__float128));

    
    for (i=0; i<klen; i++)
    {
       Pkepler.Khigh[i]=mpfr_get_float128(*(k+i),MPFR_RNDN);  
       Pkepler_high.K[i]=mpfr_get_float128(*(k+i),MPFR_RNDN);
 
    }    



/* ----------- problem parameters-----------------------------------------*/

    system.ProjFun=0;
    system.StartFun=0;
    options.TheOutput=0;
    Pkepler.RR=0;
    Pkepler.DRR=0;
    Pkepler.DK=0;
    Pkepler.KComputation=0;

    Pkepler_high.RR=0;
    Pkepler_high.DRR=0;
    Pkepler_high.DK=0;

    system.ode=ode;

    switch (ode)
    {

     case 1: 			

           system.f = NbodyOde;
           Pkepler.RR=RR1;	      
           Pkepler.DRR=DRR1; 
           Pkepler.KComputation=KComp1;     
           system.ham= HamNbody;

           system_high.f=NbodyOde_high;
           Pkepler_high.RR=RR1;	      
           Pkepler_high.DRR=DRR1_high;    
           system_high.ham= HamNbody;


     break;

     case 2: 			     
 
           system.f = NbodyOde;
           Pkepler.RR=RR2;	      
           Pkepler.DRR=DRR2; 
           Pkepler.KComputation=KComp2;     
           system.ham= HamNbody;

           system_high.f=NbodyOde_high;
           Pkepler_high.RR=RR2;	      
           Pkepler_high.DRR=DRR2_high;    
           system_high.ham= HamNbody;


     break;

     default:
           printf("error. ode\n");
     break;

     }

     
    if (mixed==true)
    {

        system.ProjFun=ProjFun_high;
        options.TheOutput=OutputFun_high;
        system.StartFun=StartFun_high;

    }
    else
    {
        system.ProjFun=ProjFun;
        options.TheOutput=OutputFun;
        system.StartFun=StartFun;

    }


    if (Pkepler.KComputation !=0)
    {
        Pkepler.KComputation(nbody,&system,&Pkepler);
        system.params.Pkepler.K=&Pkepler.K[0];
        system.params.Pkepler.K2=&Pkepler.K2[0];
        system.params.Pkepler.Khigh=&Pkepler.Khigh[0];
        system.params.Pkepler.Mu=&Pkepler.Mu[0];
        system.params.Pkepler.Muhigh=&Pkepler.Muhigh[0];

       
        system_high.params.Pkepler.K=&Pkepler_high.K[0];      
        system_high.params.Pkepler.Mu=&Pkepler.Muhigh[0];

        for (i=0; i<nbody; i++)
            Pkepler_high.K2[i]=0.;
        system_high.params.Pkepler.K2=&Pkepler_high.K2[0];    


    }
/* ----------- Integration options---------------------------------------*/


    system.params.Pkepler.RR=Pkepler.RR;
    system.params.Pkepler.DRR=Pkepler.DRR;
    system.params.Pkepler.DK=Pkepler.DK;
    system.params.Pkepler.keplerkop=Nkepler;	

    system_high.params.Pkepler.RR=Pkepler_high.RR;
    system_high.params.Pkepler.DRR=Pkepler_high.DRR;
    system_high.params.Pkepler.DK=Pkepler_high.DK;
    system_high.params.Pkepler.keplerkop=Nkepler;		
                                                       

    options.rtol=malloc(neq*sizeof(val_type));
    options.atol=malloc(neq*sizeof(val_type));

    for (i=0; i<system.neq; i++)
    {
          options.rtol[i]=RTOL;
          options.atol[i]=ATOL;
    }

    strncpy(options.filename, myfilename,STRMAX);

    options.adaptive=adaptive;


/* ----------- Integration interval----------------------------------------*/

    options.sampling=sampling;

/* ----------- Integration of the problem  ------------------------------------*/

    thread_count=threads;
    
    MallocGsCoefficients (DIR_COEFF,&gsmethod,KMAX);
    MallocGsCoefficients_high (DIR_COEFF,&gsmethod_high,KMAX);
    AssignGsCoefficients (DIR_COEFF,&gsmethod,hq,1);
    AssignGsCoefficients_high (DIR_COEFF,&gsmethod_high,hq,1);   

    method.coeffs=gsmethod;
    method.coeffs_h=gsmethod_high;

    InitStat(&system,&gsmethod,&cache_stat,&cache_vars);
    if (adaptive==true)
       InitStat_high (&system_high,&gsmethod_high,&cache_vars_high);

    ode_system.system=system;
    ode_system.system_h=system_high;

    cache.cache_stat=cache_stat;
    cache.cache_vars=cache_vars;
    cache.cache_vars_h=cache_vars_high;

    wtime0= time(NULL);
    clock0= clock();

    
    Main_FCIRK (t0q,t1q,hq, &method,&u, &ode_system, &options,&cache);

 
    clock1=clock();
    wtime1= time(NULL);


/* --------- Results ---------------------------------------------------------*/

    result_array[0]=cache.cache_stat.stepcount;
    result_array[1]=cache.cache_stat.fcn;
    result_array[2]=(float) (clock1 - clock0)/CLOCKS_PER_SEC;
    result_array[3]=(float) (wtime1 - wtime0);
    result_array[4]=cache.cache_stat.execution;
    result_array[5]=cache.cache_stat.itzero;
    result_array[6]=cache.cache_stat.totitcountzero;
    result_array[7]=cache.cache_stat.totitcount;
    result_array[8]=cache.cache_stat.MaxDE;
    result_array[9]=cache.cache_stat.hstepcount;

    free(system.params.rpar);
    free(system.params.rparhigh);
    free(system.params.ipar);

    free(u.uu);
    free(u.uul);

    free(options.rtol);
    free(options.atol);

    free(gsmethod.m);
    free(gsmethod.a);
    free(gsmethod.b);
    free(gsmethod.hb);
    free(gsmethod.c);
    free(gsmethod.hc);
    free(gsmethod.nu);
    free(gsmethod.ttau);

    free(cache_vars.z);
    free(cache_vars.li);
    free(cache_vars.fz);
    free(cache_vars.zold);
    free(cache_vars.DMin);
    free(cache_vars.delta);
    free(cache_vars.avg_delta);

    if (adaptive==true)
    {
         free(cache_vars_high.z);
         free(cache_vars_high.li);
         free(cache_vars_high.fz);
         free(cache_vars_high.zold);
         free(cache_vars_high.DMin);
    }

    free(Pkepler.K);
    free(Pkepler.K2);
    free(Pkepler.Khigh);
    free(Pkepler.Mu);
    free(Pkepler.Muhigh);

    free(Pkepler_high.K);
    free(Pkepler_high.K2);
    free(Pkepler_high.Mu);



    return 0;

}



