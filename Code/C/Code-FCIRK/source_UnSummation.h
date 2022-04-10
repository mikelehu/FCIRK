/*-----------------------------------------------------------------------------*/
/*									       */
/*                        source_UnSummation.h    	         	       */
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



   for (is=0; is<ns; is++)
    {
    
        isn=neq*is;
        
        for (i=0; i<neq; i++)
        {
            x=u->uu[i];
            xx=u->ee[i];
            y=li[isn+i];
            yy=fz[isn+i]*gsmethod->hb[is]-y;
            sub2(x,xx,y,yy,&u->uu[i],&u->ee[i]);
       }
                
    }
    


#ifdef NDEBUG

   double aux;

   aux=0.;
   for (i=0; i<neq; i++) 
   {

     aux+=u->uu[i]*u->uu[i];
   }
   printf("\nSummation Norm(w)=%lg\n",sqrt(aux));

#endif


    return;

}

