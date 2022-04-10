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
     BASE dq[dim],dv[dim];
     BASE ux[neq];

     Pkepler_sys *Pkepler;

/* ----- Implementation ------------------------------------------------------*/

     Pkepler=&params->Pkepler;

#if HIGH == 0
     val_type  x,xx,y,yy;
     k=Pkepler->K;
#else
     LOW  x,xx,y,yy;	
     k=Pkepler->K;    
#endif

     nd=neq/2;

#if HIGH == 0

     for (i=0; i<neq; i++) ux[i]=u->uu[i];

     for (i = 0; i<keplerkop; i++)
     {
           q1=i*dim;
           p1=q1+nd;

           ki=k[i];                       
           KeplerFlow(ki,&ux[q1],&ux[p1],dq,dv,h);

#else

     for (i=0; i<neq; i++)
     {      
            ux[i]=u->uu[i];
            ux[i]+=u->ee[i];
     }  

     for (i = 0; i<keplerkop; i++)
     {
           q1=i*dim;
           p1=q1+nd;

           ki=k[i];
           KeplerFlow_high(ki,&ux[q1],&ux[p1],dq,dv,h);

#endif

          for (id=0; id<dim; id++)
          {
              i1=q1+id;
              i2=p1+id;
             
              x=u->uu[i1];
              xx=u->ee[i1];
              y=dq[id];
              yy=dq[id]-y;
              add2(x,xx,y,yy,&u->uu[i1],&u->ee[i1]);
              
              x=u->uu[i2];
              xx=u->ee[i2];
              y=dv[id];
              yy=dv[id]-y;
              add2(x,xx,y,yy,&u->uu[i2],&u->ee[i2]);
              

                          
          }


     }

     return;

 }
