/*-----------------------------------------------------------------------------*/
/*									       */
/*                        source_KeplerFlowAll2.h         		       */
/*			  (stage initialization)			       */
/*									       */
/* ----------------------------------------------------------------------------*/
{

/*------ Declarations --------------------------------------------------------*/

     int dim=3;
     int i,nd,q1,p1;
     int id,i1,i2;

     BASE2B dr[dim],dv[dim];

#if HIGH2B ==0 
     Pkepler_sys *Pkepler;
#else
     Pkepler_sys_high *Pkepler;
#endif

/* ----- Implementation ------------------------------------------------------*/

     Pkepler=&params->Pkepler;


     BASE2B *k;
     k=Pkepler->K;
     nd=neq/2;

     for (i = 0; i<keplerkop; i++)
     {
          q1=i*dim;
          p1=q1+nd;

          KeplerFlow(k[i],&u[q1],&u[p1],dr,dv,h);
    
          for (id=0; id<dim; id++)
          {
              i1=q1+id;
              i2=p1+id;
              u[i1]+=dr[id];
              u[i2]+=dv[id];
          }

     }

     return;

 }
