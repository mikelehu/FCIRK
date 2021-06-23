/*------------------------------------------------------------------------------*/
/*										*/
/*         GaussCoefficients.c							*/
/*              Read in quadruple precision from a binary file,                 */						
/*		Implicit Runge-Kutta (based Gauss quadrature)  coefficients.    */
/*              s=2, s=6, s=8, s=9, s=12                                        */
/*              			                                        */
/*              	Method's coefficients: m,a,b,c.                         */
/*                      Interpolate coefficients: lambda                        */
/*                                                                              */
/*                                                                              */
/* -----------------------------------------------------------------------------*/



#include <GaussCoefficients.h>


void MallocGsCoefficients
(const char *path,tcoeffs *method, int kmax)
{

#define BASE val_type
#define HIGH 0
#include <source_MallocGsCoefficients.h>
#undef BASE
#undef HIGH

}

void MallocGsCoefficients_high
(const char *path,tcoeffs_h *method, int kmax)
{

#define BASE highprec
#define HIGH 1
#include <source_MallocGsCoefficients.h>
#undef BASE
#undef HIGH

}



void AssignGsCoefficients
(const char *path, tcoeffs *method,val_type h, int k)
{

#define BASE val_type
#define HIGH 0
#include <source_AssignGsCoefficients.h>
#undef BASE
#undef HIGH

}


void AssignGsCoefficients_high
(const char *path, tcoeffs_h *method,highprec h, int k)
{

#define BASE highprec
#define HIGH 1
#include <source_AssignGsCoefficients.h>
#undef BASE
#undef HIGH

}


void UpdateGsCoefficients_high
(tcoeffs_h *method,highprec h, int k)
{

/*------ Declarations --------------------------------------------------------*/

     int i,ii,isn;
     int ns=method->ns;



/* ----- Implementation ------------------------------------------------------*/

 

/*---- Calculate hb coefficients -----------------------------------------------*/ 
/*---- Calculate hbi coeficients in high precision ----------------------------*/

     __float128 sumq,auxq,hq;
     sumq=0.;
     hq=h;
     for (i=1; i<ns-1; i++)
     {
        auxq=method->b[i];
        method->hb[i]=hq*auxq;
        sumq+=method->hb[i];
     }

     method->hb[0]=(hq-sumq)*0.5;
     method->hb[ns-1]=(hq-sumq)*0.5;


/*---- Calculate hc coefficients -----------------------------------------------*/

     for (i=0; i<ns; i++)
        method->hc[i]=h*method->c[i];


     for (ii=0; ii<k; ii++)
     {
         isn=ii*ns;
         
         for (i=0; i<ns; i++) 
           method->ttau[isn+i]=ii*h+method->hc[i]-k*h/2;

    }

/*---- Verify that Sum hbi -h=0 ------------------------------------------------*/

#    ifdef DEBUGCOEFFICIENTS
     printf("Sum hbi -h =0");
     sumq=0.;
     hq=h;
     for (i=0; i<ns; i++)
     {
       auxq= method->hb[i];
       sumq+=auxq;
     }
     n = quadmath_snprintf(buf, sizeof buf, "%+-#*.2Qe", width, sumq-hq);
     if ((size_t) n < sizeof buf) printf("%s\n",buf);


     double aux2;
 
     for (i=0; i<ns; i++)
     {
          for (j=0; j<ns+1; j++)
          {
            aux2=method->nu[i*(ns+1)+j];
            printf("%lg,",aux2);
          }
          printf("\n");
     }
#    endif


     return;
}









