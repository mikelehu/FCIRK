/*-----------------------------------------------------------------------------*/
/*									       */
/*                        source_NbodyGFcn.h      	       		       */
/*									       */
/* ----------------------------------------------------------------------------*/
{


/*------ Declarations --------------------------------------------------------*/

#if HIGHNFCN ==0 
   Pkepler_sys *Pkepler;
#else
   Pkepler_sys_high *Pkepler;
#endif

/* ----- Implementation ------------------------------------------------------*/


     Pkepler=&params->Pkepler;
     Pkepler->DRR(neq,t,u,dR,params);

 
     return;

}
