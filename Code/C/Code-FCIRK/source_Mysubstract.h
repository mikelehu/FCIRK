/*-----------------------------------------------------------------------------*/
/*									         */
/*                        source_Mysubstract.h   	       		 */
/*									         */
/* ----------------------------------------------------------------------------*/
{
/*------ Declarations --------------------------------------------------------*/

    int dim=3;

    int id;
    BASE a3;

    BASE a,b;
    BASE l1,l2,l3,l4,l5,l6,l7,l8,l9;

/* ----- Implementation ------------------------------------------------------*/

    l1=0.;
    l2=0.;
    l7=0.;

    for (id=0; id<dim; id++)
    {
      l1+=q[id]*q[id];
      l2+=(q[id]+Dq[id])*(q[id]+Dq[id]);
      l7-=(2*q[id]+Dq[id])*Dq[id];
    }

    l3=PSQRT(l1);
    l4=PSQRT(l2);
    b=1./l3;
    a=1./l4;
    a3=1./(l4*l2);

    l5=a+b;
    l6=l5*a+1./l1;   //a^2+ab+b^2

    l8=1./(l3+l4);
    l9=l7*l8*a*b*l6;

    for (id=0; id<dim; id++) sub[id]=l9*q[id]+a3*Dq[id];


#ifdef DEBUGRR2

    printf("Mysustract\n");
    double maux;
    for (id=0; id<dim; id++)
    { 
      maux=sub[id];
      printf("%.20lg,",maux);
    
    }
    printf("\n");

#endif

    return;

}

