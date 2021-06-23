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
#if HIGH ==0
     BASE *fz;
     fz=cache_vars->fz;
#else
     BASE *li;
     li=cache_vars->li;
#endif

/* ----- Implementation ------------------------------------------------------*/


//---- High-prec computation of Li 



#if HIGH ==0
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


#ifdef NDEBUG

   double aux;
//   double aux2;

   aux=0.;
   for (i=0; i<neq; i++) 
   {
//     aux2=u->uu[i];
//     printf("%lg,",aux2);
     aux+=u->uu[i]*u->uu[i];
   }
   printf("\nSummation Norm(w)=%lg\n",sqrt(aux));

#endif


    return;

}

