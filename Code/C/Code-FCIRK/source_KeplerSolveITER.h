/*-----------------------------------------------------------------------------*/
/*									       */
/*                        source_KeplerSolveITER.h         		       */
/*									       */
/* ----------------------------------------------------------------------------*/

{
/*   G1^2+(beta *G2-2)*G2 =0,                          */
/*   G1 - Sin(x0)/nu =0,                               */
/*                                                     */
/*   gamma X + eta G2 + rkg G1 == dt -----> X, G1, G2  */

/*------ Declarations --------------------------------------------------------*/


     BASE nu,rkg;
     BASE X0,x0;
     BASE G0,G1,G2;

     val_type rh0,F1,F3,lag;
     val_type dG1,dG2,dX;


/* ----- Implementation ------------------------------------------------------*/
   

     nu=KSQRT(beta);			// high precision

     rkg=r0-gamma;			// high precision
     X0=*X;				// high precision
     x0=nu*X0;				// high precision

     G2=G[2];				// high precision
     G1=KSIN(x0)/nu;			// high precision
     G0=G[0];                           // high precision
     rh0=gamma+rkg*G0;                  // low  


     F1=(G1*G1+(beta*G2-2)*G2)/G0;      // high precision
     F3=gamma*X0+eta*G2+rkg*G1-dt;	// high precision	   
     lag=2*eta*G1+2*rh0;		// low


     dX=-(eta*F1+2*F3)/lag;		// low
     dG1=G0*dX;		         	// low
     dG2=(-2*F3*G1+F1*rh0)/lag;	        // low

     /* results */

     *X=X0+dX;				// high precision
     G[1]=G1+dG1;			// high precision
     G[2]=G2+dG2;			// high precision
     G[0]=1-beta*G[2]; 			// high precision


}
