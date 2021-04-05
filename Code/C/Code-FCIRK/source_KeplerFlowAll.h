/*-----------------------------------------------------------------------------*/
/*									       */
/*                        source_KeplerFlowAll.h         		       */
/*									       */
/* ----------------------------------------------------------------------------*/
{

/*------ Declarations --------------------------------------------------------*/

     int dim=3;
     int i,nd,q1,p1;
     int id,i1,i2;

     val_type *k;
     BASE ki;
     BASE dr[dim],dv[dim];
     BASE *ux;

     Pkepler_sys *Pkepler;


/* ----- Implementation ------------------------------------------------------*/

     Pkepler=&params->Pkepler;

#if HIGH == 0
     k=Pkepler->K;
     ux=u->uul;   //low
#else
     k=Pkepler->K;
     ux=u->uu;	//high     
#endif

     nd=neq/2;


#if HIGH == 0

     for (i = 0; i<keplerkop; i++)
     {
           q1=i*dim;
           p1=q1+nd;

           ki=k[i];

           KeplerFlow(ki,&ux[q1],&ux[p1],dr,dv,h);

#else

     for (i = 0; i<keplerkop; i++)
     {
           q1=i*dim;
           p1=q1+nd;

           ki=k[i];

           KeplerFlowMIXED_high(ki,&ux[q1],&ux[p1],dr,dv,h);
#endif

          for (id=0; id<dim; id++)
          {
              i1=q1+id;
              i2=p1+id;
              u->uu[i1]+=dr[id];
              u->uu[i2]+=dv[id];
              u->uul[i1]=u->uu[i1];
              u->uul[i2]=u->uu[i2];
          }


     }

     return;

 }
