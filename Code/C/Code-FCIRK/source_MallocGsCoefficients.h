/*-----------------------------------------------------------------------------*/
/*									         */
/*                        Julia_MallocGsCoefficients.h	       	 */
/*									         */
/* ----------------------------------------------------------------------------*/
{
/*------ Declarations --------------------------------------------------------*/

     int ns=method->ns;
     int i,info;

     char mydir[256];
     FILE *fileNU;
     char filenameNU[STRMAX];

     __float128  *nu;
     nu = (__float128 *) malloc(ns*(ns+1)*sizeof(__float128));

     method->m =(BASE *) malloc(ns*ns*sizeof(BASE));
     method->a =(BASE *) malloc(ns*ns*sizeof(BASE));
     method->b = (BASE *) malloc(ns*sizeof(BASE));
     method->hb = (BASE *) malloc(ns*sizeof(BASE));
     method->c = (BASE *) malloc((ns)*sizeof(BASE));
     method->hc = (BASE *) malloc((ns)*sizeof(BASE));
     method->nu = (BASE *) malloc(ns*(ns+1)*sizeof(BASE));
     method->ttau = (BASE *) malloc((kmax*ns)*sizeof(BASE));
     method->d = (BASE *) malloc((ns)*sizeof(BASE));

     strcpy(mydir,path);

     switch (ns)
     {

     case 8:

          strcat(mydir,"S8/");

          /*----d coefficients for adative step-------------------------------------*/
          /*--- computed by function diCoeffs!(alpha,T) implemented in Julia-------*/

          method->d[0]=-403.9205105285396;
          method->d[1]=1296.9422880742814;
          method->d[2]=-2168.250266304377;
          method->d[3]=2693.7993067876337;
          method->d[4]=-2693.7993067876614;
          method->d[5]=2168.250266304441;
          method->d[6]=-1296.9422880743396;
          method->d[7]=403.92051052856135;

     break;

     default:
          printf("Coefficients not defined\n");
     break;

     }
    
     strcpy(filenameNU,mydir);
     strcat(filenameNU,"QLamdaCoef.bin");  

     fileNU = fopen(filenameNU,"rb");
     if (fileNU == NULL) printf("Coefficients file doesnt exists\n");

     info=fread(nu, sizeof(__float128),ns*(ns+1),fileNU);
     if (info == -1) printf("Error fread command\n");
     for (i=0;i<ns*(ns+1); i++)	method->nu[i] = nu[i];     

     fclose(fileNU);
     free(nu);

     return;
}


