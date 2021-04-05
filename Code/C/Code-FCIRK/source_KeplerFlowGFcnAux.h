/*-----------------------------------------------------------------------------*/
/*									       */
/*                        source_KeplerFlowGFcnAux.h   	       		       */
/*									       */
/* ----------------------------------------------------------------------------*/
{

     int dim=3;

     int id;

     BASE r0inv,betainv;
     BASE lambda[16];
     BASE r0kg,G1b,G2b;
     BASE Rqv[2*dim],RQV[2*dim];


     if (k==0)			
     {

        for (id=0; id<dim; id++)
        {
          Gq[id]=gq[id]-dt*gv[id];		
          Gv[id]=gv[id];	
		
        }


     }
     else
     {


        for (id=0; id<dim; id++)
        {
           Rqv[id]=-gv[id];        // (Rq1,Rq2,Rq3)<---(-gv1,-gv2.-gv3)
           Rqv[dim+id]=gq[id];     // (Rv1,Rv2,Rv3)<---( gq1, gq2, gq3)

        }

    
        betainv=1/aux->beta;
        r0inv=1/aux->r0;
        G1b=0.5*(aux->X*aux->GG[0]-aux->GG[1])*betainv;    
        G2b=(0.5*aux->X*aux->GG[1]-aux->GG[2])*betainv;    

        lambda[15]=(-V[0]*Rqv[dim]-V[1]*Rqv[dim+1]-V[2]*Rqv[dim+2])*aux->rinv;
        lambda[14]=(-Q[0]*Rqv[dim]-Q[1]*Rqv[dim+1]-Q[2]*Rqv[dim+2])*aux->rinv;
        lambda[13]=(-V[0]*Rqv[0]-V[1]*Rqv[1]-V[2]*Rqv[2]);
        lambda[12]=-Q[0]*Rqv[0]-Q[1]*Rqv[1]-Q[2]*Rqv[2];
        lambda[11]=aux->b[2]*lambda[14]+aux->b[3]*lambda[15];
        lambda[10]=-aux->zeta*lambda[11]-aux->alpha*lambda[12]+aux->eta*lambda[13]-k*lambda[15];
        lambda[9]=-aux->eta*lambda[11]+aux->r0*lambda[13]-aux->alpha*lambda[14];
        r0kg=aux->r0-aux->gamma;
        lambda[8]=(lambda[9]*aux->GG[0]+lambda[10]*aux->GG[1])/(aux->gamma+r0kg*aux->GG[0]+aux->eta*aux->GG[1]);
        lambda[7]=-aux->GG[1]*lambda[8];
        lambda[6]=-aux->GG[2]*lambda[11];
        lambda[5]=(-lambda[7]-aux->X*lambda[8])*betainv;
        lambda[4]=2*(-aux->gamma*lambda[5]-aux->r0*lambda[6]+
                  (lambda[9]-r0kg*lambda[8])*G1b+(lambda[10]-aux->eta*lambda[8])*G2b);
        lambda[3]=(lambda[4]-aux->GG[2]*lambda[12]-aux->GG[1]*lambda[14])*r0inv;
        lambda[2]=(-aux->GG[1]*lambda[11]+aux->GG[2]*(lambda[13]-lambda[8]));
        lambda[1]=(-aux->alpha*lambda[3]-aux->beta*lambda[6]+lambda[7]-lambda[11]+aux->GG[1]*lambda[13])*r0inv;
 
   
        RQV[0]=-Q[0]*lambda[1]-V[0]*lambda[2]+(1+aux->b[0])*Rqv[0]+aux->b[2]*Rqv[dim];
        RQV[1]=-Q[1]*lambda[1]-V[1]*lambda[2]+(1+aux->b[0])*Rqv[1]+aux->b[2]*Rqv[dim+1];
        RQV[2]=-Q[2]*lambda[1]-V[2]*lambda[2]+(1+aux->b[0])*Rqv[2]+aux->b[2]*Rqv[dim+2];
        RQV[dim]=-Q[0]*lambda[2]+V[0]*lambda[4]+aux->b[1]*Rqv[0]+(1+aux->b[3])*Rqv[dim];
        RQV[dim+1]=-Q[1]*lambda[2]+V[1]*lambda[4]+aux->b[1]*Rqv[1]+(1+aux->b[3])*Rqv[dim+1];
        RQV[dim+2]=-Q[2]*lambda[2]+V[2]*lambda[4]+aux->b[1]*Rqv[2]+(1+aux->b[3])*Rqv[dim+2];

        for (id=0; id<dim; id++)
        {
             Gq[id]=RQV[dim+id];	 // (dQ1,dQ2,dQ3)<---(RV1,RV2.RV3)
             Gv[id]=-RQV[id];	 	 // (dV1,dV2,dV3)<---(-RQ1,-RQ2,-RQ3)
        }

#ifdef MDEBUG
        if (isnan(Gq[0]))
        {
           printf ("NAN KeplerFlowGFcnAux......\n");
           double dtx,Q1x, Q2x, Q3x;
           double V1x, V2x, V3x; 
           double gq1x,gq2x,gq3x;
           double gv1x,gv2x,gv3x;
           double kx;

           dtx=dt;
           Q1x=Q[0]; Q2x=Q[1]; Q3x=Q[2];
           V1x=V[0]; V2x=V[1]; V3x=V[2];
           gq1x=gq[0]; gq2x=gq[1]; gq3x=gq[2];
           gv1x=gv[0]; gv2x=gv[1]; gv3x=gv[2];
           kx=k;

           printf("dt=%lg, Q=%lg%lg,%lg, V=%lg,%lg,%lg\n", dtx, Q1x,Q2x,Q3x,V1x,V2x,V3x);
           printf("gq=%lg,%lg,%lg, gv=%lg,%lg,%lg, k=%lg\n",gq1x,gq2x,gq3x,gv1x,gv2x,gv3x,kx);  
        }
#endif
     }


     return;

}

