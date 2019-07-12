****************************************************************************
Readme.txt
****************************************************************************
     FCIRK implementation based on Fixed-Point iteration.

     version: 1.0 (2019-07-12).

     Article: "New Integration Methods for Perturbed ODEs Based on Symplectic
              Implicit Runge–Kutta Schemes with Application to Solar System
              Simulations" (Journal of Scientific Computing 2018). 

********************************************************************************
CONTENTS:

   def.c:               Parameters and general definitions we use in the code.
                        You must specify math-functions of your Odefun. 

   GaussCoefficients.c: Function to load Butcher tableau of the Implicit Runge-Kutta
                    method and the coefficients.We offer Gauss method coefficients for 
                        s=2, s=6, s=8, s=9, s=16.

   Julia_Source_FCIRK.h: 
        
   Common-FCIRK.c:      Numerical integration method. We will find these functions:

							
   Problems.c:  we provide ODE systems and Hamiltonians: 

        DRR1: perturbation part of N-Body problem of the Solar System
        DRR2: perturbation part of N-Body problem of the Solar System with Moon 
           as separate object.
		


*********************************************************************************
OUTPUT:

       The solution is returned in the file specified by the user.
   filename will contain (ti,ui)  values at each sampling time.


********************************************************************************


    
    



