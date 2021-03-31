/*-----------------------------------------------------------------------------*/
/*									       */
/*                        source_DRR1.h      		       		       */
/*									       */
/* ----------------------------------------------------------------------------*/
{

/*------ Declarations --------------------------------------------------------*/

     int dim,nbody;
     dim=3;
     nbody=neq/(2*dim);

     int i,id,i1,i2,j,j1,j2;
     int nd;
     BASEDRR1 d3,da,qij[dim];
     BASEDRR1 *Gm,GmSun,Gma[nbody];

     BASEDRR1 *mu,*k2;

#if HIGHDRR1 ==0 
   Pkepler_sys *Pkepler;
#else
   Pkepler_sys_high *Pkepler;
#endif

/* ----- Implementation ------------------------------------------------------*/

     nd=neq/2;
     Gm=params->rpar;
     Pkepler=&params->Pkepler;
     mu=Pkepler->Mu;
     k2=Pkepler->K2;

     GmSun=Gm[0];
     for (i=0; i<nbody; i++) {Gma[i]=Gm[i+1];}

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

     for (i=0; i<nbody; i++)
     {
          i1=i*dim;
          i2=nd+i1;

          for (j=i+1; j<nbody; j++)
          {
                j1=j*dim;
                j2=nd+j1;
                d3=0.;
                for (id=0; id<dim; id++)
                {
                     qij[id]=(u[i1+id]-u[j1+id]);
                     d3+=qij[id]*qij[id];
                }

                d3=PSQRT(d3)*d3;

                for (id=0; id<dim; id++)
                {
                     dR[i1+id]+=mu[j]*u[j2+id];
                     dR[j1+id]+=mu[i]*u[i2+id];
                     dR[i2+id]-=Gma[j]*qij[id]/d3;
                     dR[j2+id]+=Gma[i]*qij[id]/d3;
                }
           }

          for (id=0; id<dim; id++)
          {
                dR[i1+id]*=1/GmSun;
                dR[i2+id]=(Gma[i]*dR[i2+id])/mu[i];
          }

     }


     for (i=0; i<nbody; i++)
     {

         if (k2[i]!=0)
         {
              i1=i*dim;
              i2=nd+i1;

              da=0.;
              for (id=0; id<dim; id++) da+=u[i1+id]*u[i1+id];
              da=PSQRT(da)*da;
                 
              for (id=0; id<dim; id++)
                   dR[i2+id]-=k2[i]*u[i1+id]/da;
 
         }

     }
        


     return ;

}

