/*-----------------------------------------------------------------------------*/
/*									       */
/*                        source_Summation.h         	         	       */
/*			 						       */
/*									       */
/* ----------------------------------------------------------------------------*/

{


/* ----- First initializations -----------------------------------------------*/

     int neq,ns;

     neq=system->neq;
     ns=gsmethod->ns;

/*------ Declarations --------------------------------------------------------*/

     int i,is,isn;
#if HIGHSUM ==0
     BASESUM *fz;
     fz=cache_vars->fz;
#else
     BASESUM *li;
     li=cache_vars->li;
#endif

/* ----- Implementation ------------------------------------------------------*/


//---- High-prec computation of Li 


#if HIGHSUM ==0
    for (is=0; is<ns; is++)
    {
        isn=neq*is;
        u->uu[0]+=fz[isn]*gsmethod->hb[is];

        for (i = 1; i<neq; i++)
                u->uu[i]+=fz[isn+i]*gsmethod->hb[is];
    }
#else
 
   for (is=0; is<ns; is++)
    {
        isn=neq*is;
        u->uu[0]+=li[isn];

        for (i = 1; i<neq; i++)
                u->uu[i]+=li[isn+i];
    }

#endif

    for (i = 0; i<neq; i++) u->uul[i]=u->uu[i];


#ifdef MDEBUG

   double aux;

   aux=0.;
   for (i=0; i<neq; i++) aux+=u->uu[i]*u->uu[i];
   printf("Norm(w)=%lg\n",sqrt(aux));

#endif


    return;

}
