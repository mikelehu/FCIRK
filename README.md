FCIRK: Flow-Composed Implicit Runge-Kutta Method
==================================================


Flow Composed Implicit Runge-Method (FCIRK) implementation is developed in C-programming language and the integration method is interfaced via Jupyter/Julia powerful computational enviroment.

Related article: Article: "New Integration Methods for Perturbed ODEs Based on Symplectic Implicit Runge–Kutta Schemes with Application to Solar System   Simulations" (Journal of Scientific Computing 2018)".

# Repository content 

##  Documentation

*  Installation notes
*  N-Body problem equations

##  Code

      C:  C-language implementation.

            Code-Binary
            Code-FCIRK
            CoefficientsData 
	    	
       Julia (Jupyter Notebooks)

            C-Interfaces 
            Auxiliar functions
            NBodyProblem

##  Examples/16-Body (Jupyter Notebooks)
        
	1-Experiment (FCIRK-S8LQ h=3)
	2-Experiment (Efficiencuy)


#  Requirements

The current version of FCIRK software has been designed to work in Linux.

* Requires Ubuntu 12.04 or later.
* 64-bit systems have been tested successfully.
* You will need the normal gcc C/C++ compiler toolchain.

You will probably have to install a variety of packages which aren't normally installed.


# Endnotes

The project was created by Mikel Antoñana.

Date: 2019/07/12
