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
    printf("neq=%i, Nkekepler=%i, Moreq=%i\n",neq,Nkepler,Moreq);
#endif
 
/*------ declarations --------------------------------------------------------*/

    int i;

    gauss_method gsmethod;
    solution u;
    toptions options;
    ode_sys system;
    Pkepler_sys Pkepler;
    solver_stat thestat;

    clock_t clock0, clock1;
    time_t  wtime0,wtime1;

    val_type hq,t0q,t1q;

#ifdef INITVALUES
     int n;
     int width = 40;  
     char buf[128];
#endif

/* ----------- implementation  --------------------------------------------*/

    mpfr_set_default_prec (256);  // BigFloat's default precision in Julia


/* ----------- Initial values--------------------------------------------*/

    system.neq=neq;
    gsmethod.ns=ns;

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
    system.params.ipar =(int *)malloc(ilen*sizeof(int));
    system.params.numrpar=rlen;
    system.params.numipar=ilen;

    for (i=0; i<rlen; i++)
    {
         system.params.rparhigh[i]=mpfr_get_float128(*(rpar+i),MPFR_RNDN);
         system.params.rpar[i]=system.params.rparhigh[i];
    }


#ifdef INITVALUES
     printf("Rpar\n");
     for(i=0; i<rlen; i++)
     {
       n = quadmath_snprintf(buf, sizeof buf, "%+-#*.36Qe", width, system.params.rparhigh[i]);
       if ((size_t) n < sizeof buf) printf("%s\n",buf);
     }
#endif

 
    for (i=0; i<ilen; i++) 
          system.params.ipar[i]=ipar[i];           

    Pkepler.K =(val_type *)malloc(nbody*sizeof(val_type));
    Pkepler.K2 =(val_type *)malloc(nbody*sizeof(val_type));
    Pkepler.Khigh =(__float128 *)malloc(nbody*sizeof(__float128));
    Pkepler.Mu =(val_type *)malloc(nbody*sizeof(val_type));
    Pkepler.Muhigh =(__float128 *)malloc(nbody*sizeof(__float128));
    
    for (i=0; i<klen; i++)
       Pkepler.Khigh[i]=mpfr_get_float128(*(k+i),MPFR_RNDN);      



/* ----------- problem parameters-----------------------------------------*/

    system.ProjFun=0;
    system.StartFun=0;
    options.TheOutput=0;
    Pkepler.RR=0;
    Pkepler.DRR=0;
    Pkepler.DK=0;
    Pkepler.KComputation=0;

    system.codfun=codfun;

    switch (codfun)
    {


/*
     Perturbation DRR1, valid for:
            6-body  problem
            10-body problem (Earth+Moon one body)
*/

     case 02: 			       // FCIRK 

           system.ProjFun=ProjFun;
           options.TheOutput=OutputFun;
           system.StartFun=StartFun;

           system.f = NbodyOde;
           Pkepler.RR=RR1;	      
           Pkepler.DRR=DRR1; 
           Pkepler.KComputation=KComp1;     
           system.ham= Ham3;

     break;


     case 03: 			       // FCIRK (Mixed-Precision)

           system.ProjFun=ProjFun_high;
           options.TheOutput=OutputFun_high;
           system.StartFun=StartFun_high;

           system.f = NbodyOde;
           Pkepler.RR=RR1;	          
           Pkepler.DRR=DRR1; 
           Pkepler.KComputation=KComp1;                
           system.ham= Ham3;

     break;


/*
     Perturbation DRR2, valid for (with Moon as separate body):
           11-body problem
           16-body problem

*/


     case 12: 			      // FCIRK 
 
           system.ProjFun=ProjFun;
           options.TheOutput=OutputFun;
           system.StartFun=StartFun;

           system.f = NbodyOde;
           Pkepler.RR=RR2;	      
           Pkepler.DRR=DRR2; 
           Pkepler.KComputation=KComp2;     
           system.ham= Ham3;

     break;




     case 13: 			       // FCIRK (Mixed-Precision)

           system.ProjFun=ProjFun_high;
           options.TheOutput=OutputFun_high;
           system.StartFun=StartFun_high;

           system.f = NbodyOde;
           Pkepler.RR=RR2;	          
           Pkepler.DRR=DRR2; 
           Pkepler.KComputation=KComp2;                
           system.ham= Ham3;

     break;


     default:
           printf("error. codfun\n");
     break;

     }

    if (Pkepler.KComputation !=0)
    {
        Pkepler.KComputation(nbody,&system,&Pkepler);
        system.params.Pkepler.K=&Pkepler.K[0];
        system.params.Pkepler.K2=&Pkepler.K2[0];
        system.params.Pkepler.Khigh=&Pkepler.Khigh[0];
        system.params.Pkepler.Mu=&Pkepler.Mu[0];
        system.params.Pkepler.Muhigh=&Pkepler.Muhigh[0];

    }
/* ----------- Integration options---------------------------------------*/


    system.params.Pkepler.RR=Pkepler.RR;
    system.params.Pkepler.DRR=Pkepler.DRR;
    system.params.Pkepler.DK=Pkepler.DK;
    system.params.Pkepler.keplerkop=Nkepler;		
                                                       

    options.rtol=malloc(neq*sizeof(val_type));
    options.atol=malloc(neq*sizeof(val_type));

    for (i=0; i<system.neq; i++)
    {
          options.rtol[i]=RTOL;
          options.atol[i]=ATOL;
    }

    strncpy(options.filename, myfilename,STRMAX);


/* ----------- Integration interval----------------------------------------*/

    options.sampling=sampling;

/* ----------- Integration of the problem  ------------------------------------*/

    thread_count=threads;
    GaussCoefficients(DIR_COEFF,&gsmethod,hq);    
    InitStat(&system,&gsmethod,&thestat);

    wtime0= time(NULL);
    clock0= clock();
    
    Main_FCIRK (t0q,t1q,hq, &gsmethod, &u, &system, &options, &thestat);
 
    clock1=clock();
    wtime1= time(NULL);


/* --------- Results ---------------------------------------------------------*/

    result_array[0]=thestat.stepcount;
    result_array[1]=thestat.fcn;
    result_array[2]=(float) (clock1 - clock0)/CLOCKS_PER_SEC;
    result_array[3]=(float) (wtime1 - wtime0);
    result_array[4]=thestat.execution;
    result_array[5]=thestat.itzero;
    result_array[6]=thestat.totitcountzero;
    result_array[7]=thestat.totitcount;
    result_array[8]=thestat.MaxDE;

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

    free(thestat.z);
    free(thestat.li);
    free(thestat.fz);

    free(thestat.zold);

    free(Pkepler.K);
    free(Pkepler.K2);
    free(Pkepler.Khigh);
    free(Pkepler.Mu);
    free(Pkepler.Muhigh);


    return 0;

}



