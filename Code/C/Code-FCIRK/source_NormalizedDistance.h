/*-----------------------------------------------------------------------------*/
/*									       */
/*                        source_NormalizedDistance.h         		       */
/*			 						       */
/*									       */
/* ----------------------------------------------------------------------------*/

{
/*---------------- Declarations ----------------------------------------------*/

     int i,is,ix;
     val_type maxi,mayi,relerrors;
     val_type maxz,maxzold;

/* --------------- Implementation --------------------------------------------*/

     maxi=0.;

     for (i=0; i<neq; i++)
     {
          maxz=0.;
          maxzold=0.;

          for (is=0; is<ns;is++)
          {
               ix=neq*is+i;
               maxz=FMAX(maxz,FABS(z[ix]));
               maxzold=FMAX(maxzold,FABS(zold[ix]));
          }

          mayi=(maxz+maxzold)/2;
          relerrors=0.;

          for (is=0; is<ns; is++)
          {
               ix=neq*is+i;
               relerrors=FMAX(relerrors,
                              FABS(z[ix]-zold[ix])/
                                   (mayi*options->rtol[i]+options->atol[i]));
          }

          maxi=FMAX(maxi,relerrors);

     }

     return maxi;

}




