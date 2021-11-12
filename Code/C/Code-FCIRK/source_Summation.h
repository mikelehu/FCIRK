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
     
     BASE *fz,*li;
     fz=cache_vars->fz;
     li=cache_vars->li;

/*------ Declarations --------------------------------------------------------*/

     int i,is,isn;
#if HIGH ==0
     val_type  x,xx,y,yy;

#else
     LOW  x,xx,y,yy;	
#endif

/* ----- Implementation ------------------------------------------------------*/


//#if HIGH ==0
//    for (is=0; is<ns; is++)
//    {
    
//        isn=neq*is;
//        u->uu[0]+=fz[isn]*gsmethod->hb[is];

//        for (i = 1; i<neq; i++)
//                u->uu[i]+=fz[isn+i]*gsmethod->hb[is];
                              
//    }
//#else
 
//   for (is=0; is<ns; is++)
//    {    
//        isn=neq*is;
        
//        u->uu[0]+=li[isn];
//        for (i = 1; i<neq; i++)
//             u->uu[i]+=li[isn+i];                               
//    }

//#endif

   for (is=0; is<ns; is++)
    {
    
        isn=neq*is;
        
        for (i=0; i<neq; i++)
        {
            x=u->uul[i];
            xx=u->ee[i];
            y=li[isn+i];
            yy=fz[isn+i]*gsmethod->hb[is]-y;
            add2(x,xx,y,yy,&u->uul[i],&u->ee[i]);
       }
                
    }
    

//    for (i = 0; i<neq; i++) u->uul[i]=u->uu[i];
//    for (i = 0; i<neq; i++) u->uu[i]=u->uul[i];    // behin behineko 2021-11-03

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

