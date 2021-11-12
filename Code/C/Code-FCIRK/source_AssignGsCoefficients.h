/*-----------------------------------------------------------------------------*/
/*									       */
/*                        source_AssignGsCoefficients.h	       		       */
/*									       */
/* ----------------------------------------------------------------------------*/
{
/*------ Declarations --------------------------------------------------------*/

     int i,j,ii,isn;
     int info;
     int ns=method->ns;

     char mydir[256];
     FILE *fileM,*fileA,*fileB,*fileC; 

     char filenameM[STRMAX];
     char filenameA[STRMAX];
     char filenameB[STRMAX];
     char filenameC[STRMAX];

#    ifdef XDEBUG
     double aux;
#    endif

     __float128 *m, *a, *b, *c; // *nu;

     m =(__float128 *) malloc(ns*ns*sizeof(__float128));
     a =(__float128 *) malloc(ns*ns*sizeof(__float128));
     b = (__float128 *) malloc(ns*sizeof(__float128));
     c = (__float128 *) malloc((ns)*sizeof(__float128));


/* ----- Implementation ------------------------------------------------------*/

     strcpy(mydir,path);

     switch (ns)
     {
     case 6:

           strcat(mydir,"S6/");

     break;

     case 8:

           strcat(mydir,"S8/");


     break;

     case 9:

           strcat(mydir,"S9/");

     break;

     case 16:

           strcat(mydir,"S16/");

     break;

     case 2:

           strcat(mydir,"S2/");

     break;


     default:
          printf("Coefficients not defined\n");
     break;

     }

     strcpy(filenameM,mydir);
     strcpy(filenameA,mydir);
     strcpy(filenameB,mydir);
     strcpy(filenameC,mydir);


     strcat(filenameM,"QMCoef.bin");
     strcat(filenameA,"QACoef.bin");
     strcat(filenameB,"QBCoef.bin");
     strcat(filenameC,"QCCoef.bin");  


    fileM = fopen(filenameM,"rb");
    if (fileM == NULL) printf("Coefficients file doesnt exists\n");
    fileA = fopen(filenameA,"rb");
    if (fileA == NULL) printf("Coefficients file doesnt exists\n");
    fileB = fopen(filenameB,"rb");
    if (fileB == NULL) printf("Coefficients file doesnt exists\n");
    fileC = fopen(filenameC,"rb");
    if (fileC == NULL) printf("Coefficients file doesnt exists\n");


    info=fread(m, sizeof(__float128),ns*ns,fileM);
    if (info == -1) printf("Error fread command\n");
    info=fread(a, sizeof(__float128),ns*ns,fileA);
    if (info == -1) printf("Error fread command\n");
    info=fread(b, sizeof(__float128),ns,fileB);
    if (info == -1) printf("Error fread command\n");
    info=fread(c, sizeof(__float128),ns,fileC);
    if (info == -1) printf("Error fread command\n");

 
    for (i=0; i<ns; i++) method->m[i*ns+i]=m[i*ns+i];
    for (i=1; i<ns; i++)
    for (j=0; j<i; j++) method->m[i*ns+j]=m[i*ns+j];
    for (i=1; i<ns; i++)
    for (j=0; j<i; j++) method->m[j*ns+i]=1-method->m[i*ns+j];


    for (i=0;i<ns*ns; i++) method->a[i] = a[i];


    for (i=0;i<ns; i++)
    {
	method->c[i]= c[i];
	method->b[i] = b[i];
    }


/*---- Verify symplectic condition ---------------------------------------------*/

#    ifdef XDEBUG
     for (i=0; i<ns; i++)
     {
        printf("\n");
        
        for (j=0; j<ns; j++) 
        {
             aux=method->m[i*ns+j]+method->m[j*ns+i]-1.;
             printf ("%lg,", aux);
        }
     }

     printf("\n");
#    endif

/*---- Calculate hb coefficients -----------------------------------------------*/ 
/*---- Calculate hbi coeficients in high precision ----------------------------*/

     __float128 sumq,auxq,hq;
     sumq=0.;
     hq=h;
     for (i=1; i<ns-1; i++)
     {
        auxq=b[i];
        method->hb[i]=hq*auxq;
        sumq+=method->hb[i];
     }

     method->hb[0]=(hq-sumq)*0.5;
     method->hb[ns-1]=(hq-sumq)*0.5;


/*---- Calculate hc coefficients -----------------------------------------------*/

     for (i=0; i<ns; i++)
        method->hc[i]=h*method->c[i];



     for (ii=0; ii<k; ii++)
     {
         isn=ii*ns;
         
         for (i=0; i<ns; i++) 
           method->ttau[isn+i]=ii*h+method->hc[i]-k*h/2;

    }

/*---- Verify that Sum hbi -h=0 ------------------------------------------------*/

#    ifdef DEBUGCOEFFICIENTS
     printf("Sum hbi -h =0");
     sumq=0.;
     hq=h;
     for (i=0; i<ns; i++)
     {
       auxq= method->hb[i];
       sumq+=auxq;
     }
     n = quadmath_snprintf(buf, sizeof buf, "%+-#*.2Qe", width, sumq-hq);
     if ((size_t) n < sizeof buf) printf("%s\n",buf);


     double aux2;
 
     for (i=0; i<ns; i++)
     {
          for (j=0; j<ns+1; j++)
          {
            aux2=method->nu[i*(ns+1)+j];
            printf("%lg,",aux2);
          }
          printf("\n");
     }
#    endif


     fclose(fileM);
     fclose(fileA);
     fclose(fileB);
     fclose(fileC);

     free(m);
     free(a);
     free(b);
     free(c);

     return;
}



