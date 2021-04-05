/*-----------------------------------------------------------------------------*/
/*									       */
/*                        source_DRR2.h      		       		       */
/*									       */
/* ----------------------------------------------------------------------------*/
{

/* ---------- First initializations ------------------------------------*/

     int dim,nbody;
     dim=3;
     nbody=neq/(2*dim);

/*------ declarations -------------------------------------------------*/

     int i,id,i1,i2,j,j1,j2;
     int iM,iM1,iM2,iE,iE1,iE2;
     int nd;
     BASE aux1,aux2;
     BASE qa[dim],qb[dim];
     BASE qEM[dim],DqEM[dim];

     BASE da,db;
     BASE *Gm,GmSun,Gma[nbody];

     BASE GQN[dim];
     BASE sub1[dim],sub2[dim],sub3[dim];
     BASE *mu,*k2;

#if HIGH ==0 
     Pkepler_sys *Pkepler;
#else
     Pkepler_sys_high *Pkepler;
#endif

/* ----------- implementation  ---------------------------------------*/

     nd=neq/2;
     Gm=params->rpar;
     Pkepler=&params->Pkepler;
     mu=Pkepler->Mu;
     k2=Pkepler->K2;

     iM=(nbody-1);  	       // Moon
     iM1=iM*dim;
     iM2=nd+iM1;
     iE=(nbody-2);            // Earth
     iE1=iE*dim;
     iE2=nd+iE1;


     GmSun=Gm[0];


     for (i=0; i<nbody; i++) {Gma[i]=Gm[i+1];}
     for (id=0; id<dim; id++) GQN[id]=-(u[iM1+id]/Gma[iE])*Gma[iM];


     for (id=0; id<dim; id++)
     {
          DqEM[id]=(u[iM1+id]/Gma[iE])*(Gma[iM]+Gma[iE]);
          qEM[id]=-(u[iE1+id]+u[iM1+id]);
     }



     for (i=0; i<nbody; i++)
     {
          i1=i*dim;
          i2=nd+i1;
          for (id=0; id<dim; id++)
          {
               dR[i1+id]=0.;
               dR[i2+id]=0.;
          }
     }



     /* First DQi  */


     for (i=0; i<(nbody-1); i++)
     {
          i1=i*dim;
          i2=nd+i1;

          for (j=0; j<(nbody-1); j++)
          {
               j1=j*dim;
               j2=nd+j1;
  
               if (j!=i)
                   for (id=0; id<dim; id++)
                        dR[i1+id]+=(u[j2+id]*mu[j]);
           }

          for (id=0; id<dim; id++)
                 dR[i1+id]*=1/GmSun;

     }



     /* Second DVi  */



     for (i=0; i<(nbody-2); i++)
     {
           i1=i*dim;
           i2=nd+i1;

           for (j=0; j<(nbody-2); j++)
           {
                j1=j*dim;
                j2=nd+j1;

                if (j!=i)
                {
                    da=0;
                    for (id=0; id<dim; id++)
                    {
                         qa[id]=(u[i1+id]-u[j1+id]);
                         da+=qa[id]*qa[id];
                    }

                    da=PSQRT(da)*da;

                    for (id=0; id<dim; id++)
                         dR[i2+id]-=Gma[j]*qa[id]/da;
                }
           }


           da=0.;
           db=0.;
           for (id=0; id<dim; id++)
           {
                 qa[id]=u[i1+id]-u[iE1+id]-u[iM1+id];
                 qb[id]=u[i1+id]-u[iE1+id]-GQN[id];     
                 da+=qa[id]*qa[id];
                 db+=qb[id]*qb[id];
           }

           da=PSQRT(da)*da;
           db=PSQRT(db)*db;

           for (id=0; id<dim; id++)
           {
                 aux1=Gma[iM]/da*qa[id]+Gma[iE]/db*qb[id];
                 aux2=Gma[iM]*(qa[id]/da-qb[id]/db);
                 dR[i2+id]-=aux1;
                 dR[iE2+id]+=Gma[i]*aux1;
                 dR[iM2+id]+=Gma[i]*aux2;
           }

           for (id=0; id<dim; id++)
                 dR[i2+id]=(dR[i2+id]/mu[i])*Gma[i];
     }


     /* Dvi i=N-1, i=N */


#if HIGH ==0 
     Mysubstract(&u[iE1],&u[iM1],sub1);
     Mysubstract(&u[iE1],GQN,sub2);
     Mysubstract(qEM,DqEM,sub3);
#else
     Mysubstract_high(&u[iE1],&u[iM1],sub1);
     Mysubstract_high(&u[iE1],GQN,sub2);
     Mysubstract_high(qEM,DqEM,sub3);
#endif




     for (id=0; id<dim; id++)
     {
          dR[iE2+id]-=GmSun*(Gma[iE]*sub2[id]+Gma[iM]*sub1[id]);
          dR[iM2+id]+=-GmSun*(Gma[iM]*sub3[id]);
          dR[iE2+id]*=1/mu[iE];
          dR[iM2+id]*=1/mu[iM];

     }



     for (i=0; i<nbody; i++)
     {
          i1=i*dim;
          i2=nd+i1;
   
          da=0.;
          for (id=0; id<dim; id++) da+=u[i1+id]*u[i1+id];
          da=PSQRT(da)*da;
                 
          for (id=0; id<dim; id++)
               dR[i2+id]-=k2[i]*u[i1+id]/da;

     }



#ifdef MDEBUG

    double myaux;
    printf("\nDRR2B\n");
    printf("q :");    
    for (i=0; i<neq/2;i++)
    {
         myaux=dR[i];
         printf("%lg,", myaux);
    }
    printf("\n");
    printf("v :");
    for (i=neq/2; i<neq;i++) 
    {
        myaux=dR[i]; 
        printf("%lg,", myaux);
    }
    printf("\n");


#endif

     return ;

}



