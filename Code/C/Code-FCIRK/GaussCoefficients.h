/*-----------------------------------------------------------------------------*/
/*									       */
/*                        GaussCoefficients.h	         		       */
/*									       */
/* ----------------------------------------------------------------------------*/


#include <def.h>
#include <math.h>
#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <quadmath.h>



void MallocGsCoefficients
(const char *path,tcoeffs *method, int kmax
);

void MallocGsCoefficients_high
(const char *path,tcoeffs_h *method, int kmax
);

void AssignGsCoefficients
(const char *path,tcoeffs *method,val_type h, int k
);

void AssignGsCoefficients_high
(const char *path, tcoeffs_h *method,highprec h, int k
);

void UpdateGsCoefficients_high
(tcoeffs_h *method,highprec h, int k
);
