****************************************************************************
Readme.txt
****************************************************************************
     Gauss Implicit Runge-Kutta implementation. 
     Fixed-Point iteration.

     version: 1.1 (21-04-2017).

     Article: "Reducing and monitoring round-off error propagation for sympletic
               implicit Runge-Kutta schemes" (2017). 

********************************************************************************
CONTENTS:

   def.c:               Parameters and general definitions we use in the code.
                        You must specify math-functions of your Odefun. 

   Terminal-IRKFP.c:    An example to show how to run an integration from the Unix terminal.
   math-IRKFP.c:	An example to show how to run an integration from mathematica. 

   GaussCoefficients.c: Function to load Butcher tableau of the Implicit Runge-Kutta method and 
                        the coefficients of the  reformulation as explained in the article.
                        We offer Gauss method coefficients for s=6, s=8, s=16.

   Common-IRKFP.c:      Numerical integration method. We will find these functions:

	IRKFP () : 	         IRK Fixed-Point integration (main function).
	Fixed_point_Step ():     Fixed iteration Step.	
	General_FP_It():         General fixed point iteration.	
	Partitioned_FP_It():     Partitioned fixed point iteration.		
        Default_Stage_init:      Initialization of Stages, Y_{n,i}=y_n	
        Interpolated_Stage_init	 Initialization of stages, Y_{n,i}=G(Y_{n-1,i})		
											
	MyOutput():	         An output function defined by the user.

	NormalizedDistance():    Check convergence of the iteration.
	StopCriterion ():        Stopping criterion.
	RemoveDigits():          Rounding a floating point number with p-r significant binary digits.
        CompensatedSummation():  Function to compute y_{n+1}=y_{n}+ sum_{i} L_{n,i}				
							
   Problems.c:  Examples of ODE systems and Hamiltonians: Double pendulum and N-Body problems.
														
	Ode1()=OdePendulum():	
	Ham1()=HamPendulum():	
	Ode2()=OdeNBody():
	Ham2()=HamNBody():

   math-IRKFP.tm:		 Mathlink file.
			

*********************************************************************************
USE OF THE CODE:

   The Function call to integrate an ODE system is:

     IRKFP (t0,t1,h,&method,&u,&system,&options,&thestat);
      
       t0,t1:     interval of numerical integration.
       h:         step size.
       method:    structure with information about the method to be used in the integration: 
                  Butcher tableau and the reformulation of the symplectic IRK scheme (see below)
       u:         initial values of the problem (structure with two parts: uu for initial values 
                  and ee for accumulated round-off error values) 
       system:    structure with the information of the problem to be solved (function that evaluates 
                  the differential equations, the Hamiltonian, list of parameters of the problem...)
       options:   options of the integration: tolerances, function to be called after each step, 
                  name of the file where save intermediate values, initialization of the stages at each step...
       thestat:   internal data of the integration: number of iterations, number of steps, number of 
                  calls to the ODE function...

    All these structures are defined in the file "def.h". Let's see them:
       method:    It defines the scheme used to solve the equations. Its members are:
                      int ns;        the number of stages 
                      val_type *c,*b,*a;	The Butcher tableau of the method.
                      val_type *m;   mij=aij/bj coefficients derived from the reformulation to maintain symplecticity mij+mji-1=0.
                      val_type *hc;  hc=h*c (to avoid recalculations).
                      val_type *hb;  hb=h*b (to avoid recalculations).
                      val_type *nu;  coefficients for interpolation (a*/bj).
                      int *orderedindices;   ascending order indices for bi coefficients                  
                  The function defined in the file GaussCoefficients.c loads all these values for 
                  Gauss collocation schemes for s=6, s=8 and s=16. The informations needed by this function 
                  are the number of stages and the step size h.
       system:    This structure has the information of the ODE: 
                      int neq:       Number of equations or dimension of the problem. 
                      void (*f)():   the function that evaluates the differential equations, The user has to code it and set here its name.
                      val_type (*Ham)(): the function that evaluates the Hamiltonian (if the problem is Hamiltonian), 
                      parameters:    structure with the parameters that define the problem (list of real parameters and list of integer parameters),
                      int cod[2]:    if the problem is a partitioned problem there is the possibility to define which part must be evaluated first.
       options:   The user can set several options for the integration process:
                      val_type *rtol,*atol;    relative and absolute tolerances	
                      void (*TheOutput)();     The function that will be called after each step (or after "sampling" steps).	  
                      int sampling;            For long time integrations TheOutput function will be called once after "sampling" steps 
                      int rdigits,mrdigits;    When we want a computation with r digits less of accuracy we can set rdigits = r
                      int (*iteration)();      The function that will be called for fixed point iteration. 
                                               We offer one special function for partitioned problems (Partitioned_FP_It) and another function for general problems (General_FP_It).
                      char filename[STRMAX];   Output filename.
                      void (*StageInitFn)();   The function that initializes the stages at the beginning of the step. We offer two functions, but the user can set its own function. 

*********************************************************************************
OUTPUT:

       the value at t1 is returned in the variable u.
       If the user uses MyOutput function then the file set in options.filename will contain the values at each sampling time.


*********************************************************************************
PARAMETERS (file: def.h):

   You can specify next parameters:

   PARALLEL: define it to run a parallel execution. In this case we have to set the number of 
             threads in the global variable named thread_count:

             int thread_count = 1;   //  number of threads for parallel computations.

   MAXIT:  maximum number of fixed point iterations 
   RTOL,ATOL: fixed point iteration default tolerances. 


   #define DIR_TERM :  // Path Coefficients for terminal executions  
   #define DIR_MATH :  // Path Coefficients for mathematica executions



********************************************************************************
INSTALATION (Ubuntu 16.04 lts):

   Two ways to execute :

     Unix terminal execution:
          make term-IRKFP
          ./Terminal-IRKFP.exe

      Mathematica (double)
          make math-IRKFP
          Execution from mathematica notebook (see examples).


********************************************************************************
EXAMPLES:

    NCDP: Mathematica notebook integration of Double Pendulum problem (non-chaotic).
    CDP : Mathematica notebook integration of Double Pendulum problem (chaotic).
    OSS : Mathematica notebook integration of Outer solar system problem.

*********************************************************************************

    
    



