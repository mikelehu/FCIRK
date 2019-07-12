
/*****************************************************************************/
/*  File: CBinFilesWRA.c						     */
/*  Content:					   			     */
/*        CReadBin : read solution from a binary file			     */
/*        CWriteBin:  write solution in a binary file			     */
/*        CAppendBin: append solution in a binary file			     */
/*									     */
/*  File format:                                                             */
/*                                                                           */
/*       t_0,u_0                                                             */
/*       t_1,u_1                                                             */
/*       ...                                                                 */
/*       t_n,u_n                                                             */
/*                                                                           */
/*       where                                                               */
/*	   t_i: real value						     */
/*	   u_i:	vector of size neq.					     */
/*									     */
/*       All data are given in quadruple precision                           */
/*									     */
/*****************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include <string.h>
#include <quadmath.h>
#include <gmp.h>
#include <mpfr.h>

#define STRMAX 256		// Filename maximum string.


/*----------------------------------------------------------------------------*/
/*									      */
/*  CReadBin								      */
/*									      */	
/*----------------------------------------------------------------------------*/
int CReadBin (int nout, int neq, char *myfilename, mpfr_ptr *result_array)

{

/*------ ---------------------------------------------------------------------*/


    mpfr_set_default_prec (256);  

/*------ declarations --------------------------------------------------------*/

    int i,in,info;

    int k;
    FILE *ptr_myfile;
    char filename[STRMAX]; 


    /* Our structure */

    struct rec
    {
      __float128 t;
      __float128 uu[neq];
    };

    struct rec my_record;


/* ----------- implementation  -----------------------------------------------*/

    strncpy(filename,myfilename,STRMAX);
    ptr_myfile=fopen(filename,"rb");

    if (!ptr_myfile)
    {
	printf("Unable to open file!");
	return 1;
    }

    for ( k=0; k < nout; k++)
    {
	info=fread(&my_record,sizeof(struct rec),1,ptr_myfile);
        if (info == -1) printf("Error fread command\n");

        in=k*(neq+1);
        mpfr_set_float128 (*(result_array+in),my_record.t,MPFR_RNDN);

        for (i=0; i<neq; i++)
             mpfr_set_float128 (*(result_array+(in+i+1)),my_record.uu[i],MPFR_RNDN);

#ifdef DEBUG
                double aux;
                aux=my_record.t;
                printf("%lg\n",aux);
                
#endif
    }

    fclose(ptr_myfile);
    return 0;

}

/*----------------------------------------------------------------------------*/
/*									      */
/*  CwriteBin								      */
/*									      */	
/*----------------------------------------------------------------------------*/

int CWriteBin (int nout, int neq, char *myfilename, mpfr_ptr *y0)

{

/*------ ---------------------------------------------------------------------*/

   mpfr_set_default_prec (256);  

/*------ declarations --------------------------------------------------------*/

    int i,in;

    int k;
    FILE *ptr_myfile;
    char filename[STRMAX]; 


    /* Our structure */

    struct rec
    {
      __float128 t;
      __float128 uu[neq];
    };

    struct rec my_record;


/* ----------- implementation  -----------------------------------------------*/

    strncpy(filename,myfilename,STRMAX);
    ptr_myfile=fopen(filename,"wb");

    if (!ptr_myfile)
    {
	printf("Unable to open file!");
	return 1;
    }

    for ( k=0; k < nout; k++)
    {

         in=k*(neq+1);
         my_record.t=mpfr_get_float128(*(y0+in),MPFR_RNDN);

         for (i=0; i<neq; i++)
               my_record.uu[i]=mpfr_get_float128(*(y0+in+i+1),MPFR_RNDN); 


         fwrite(&my_record, sizeof(struct rec), 1, ptr_myfile);

    }


    fclose(ptr_myfile);
    return 0;

}


/*----------------------------------------------------------------------------*/
/*									      */
/*  CAppendBin								      */
/*									      */	
/*----------------------------------------------------------------------------*/
int CAppendBin (int nout, int neq, char *myfilename, mpfr_ptr *y0)

{

/*------ ---------------------------------------------------------------------*/

    mpfr_set_default_prec (256);  

/*------ declarations --------------------------------------------------------*/

    int i,in;

    int k;
    FILE *ptr_myfile;
    char filename[STRMAX]; 


    /* Our structure */

    struct rec
    {
      __float128 t;
      __float128 uu[neq];
    };

    struct rec my_record;


/* ----------- implementation  -----------------------------------------------*/

    strncpy(filename,myfilename,STRMAX);
    ptr_myfile=fopen(filename,"ab");

    if (!ptr_myfile)
    {
	printf("Unable to open file!");
	return 1;
    }

    for ( k=0; k < nout; k++)
    {
          in=k*(neq+1);
          my_record.t=mpfr_get_float128(*(y0+in),MPFR_RNDN);

          for (i=0; i<neq; i++)
               my_record.uu[i]=mpfr_get_float128(*(y0+in+i+1),MPFR_RNDN); 

          fwrite(&my_record, sizeof(struct rec), 1, ptr_myfile);

    }


    fclose(ptr_myfile);
    return 0;

}
