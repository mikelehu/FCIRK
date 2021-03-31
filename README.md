# FCIRK16
- Flow-Composed Implicit Runge-Kutta 16th order
(implementation in **C language**)
- Experiments are provided in **Jupyter Notebooks** interfaced via **Julia**



## Description

We present a **C** implementation of a 16th order Flow-Composed Implicit Runge-Kutta integrator
FCIRK16 appropiate for the numerical integration of the Solar System, or other systems with
near-Keplerian motions.

<img src="https://latex.codecogs.com/svg.latex?\&space;
\frac{du}{dt}=k(u)+g(u), \qquad u(t_0)=u_0 \in \mathbb{R}^D" />


where <img src="https://latex.codecogs.com/svg.latex?\&space;g(u) \quad "/> is considered as a perturbation, and the unperturbed system

<img src="https://latex.codecogs.com/svg.latex?\&space;\frac{du}{dt}=k(u) \quad"  /> can be solved exactly for any initial value.


## Installation

We do all of our development in Ubuntu.
You can install FCIRK16 application  but previously, you must install some C specific libraries.

### C  libraries (prerequisites)

We use some specifics libraries that you must be installed:

(1) uuid-dev

![Install uuid](/install_uuid.png)

(2) quadmath library

Quad-Precision Math Library Application Programming Interface (API).
https://gcc.gnu.org/onlinedocs/libquadmath/#toc-Typedef-and-constants-1

(3) mpfr library

MPFR is a portable library written in C for arbitrary precision arithmetic on floating-point numbers.
The MPFR library is already installed on some GNU/Linux distributions and “How to Install”
instructions, are explained in “GNU MPFR” manual.

(4) OpenMP Application Programming Interface

OpenMP API for parallelism in C, C++ and Fortran programs. Most of Compilers support Open
API. http://www.openmp.org


### FCIRK: compile and link



- libFCIRKDOUBLE.so, libFCIRKLONG.so and libFCIRKQuad.so

```
$ make -s -C ./Code/C/Code-Binary/
```
- libCBinary.so: binary format files read/write functions

```
$ make -s -C ./Code/C/Code-FCIRK/
```



### Jupyter/Julia computational enviroment

(1) Julia: Julia 1.5.3 version (Nov 9, 2020) (http://julialang.org).

(2) Jupyter:
we build our experiments on Jupyter notebooks (an open source tool for interative computing)
and calling to FCIRK solver using Julia programming language

  https://github.com/JuliaLang/IJulia.jl


## FCIRK by Example

We consider a Newtonian point-mass 16-body model of the Solar System applied
with **Heliocentric Canonical Coordinates** for all bodies except **Geocentric Coordinates** for the Moon. It includes:
- the Sun, the eight planets, Pluto,
- the Moon as a separate body, and
- the five main bodies of the asteroid belt (Ceres, Pallas, Vesta, Iris and Bamberga),

We consider the initial values at Julian day (TDB) 2440400.5 (the 28th of June of 1969) from DE430 Ephemerides,
renormalized so that the center of mass of the 16 bodies is at rest

### Step 1: Defining  the problem

To solve this numerically, we define a problem type by giving it the equation of perturbation part (see file: Problems.c):


```C
void DRR2 (int neq, val_type t,val_type *u,val_type *dR,parameters *params)
{
     < CODE >
     return ;

}
```

Define the problem giving  initial
conditions and the timespan to solve over

```julia
(u0,k,Gm)= Initial_Values_N16()
U0=ChangeBartoHel(u0,Gm)
h=3.
t0=0.
tend=100*h
tspan=(t0,tend)
ode=2
prob=ODEProblem(ode,tspan,U0,k)
```

### Step 2: Solving the problem


```julia
sol1="./sol1.bin"
destats=solve(sol1,prob,longFCIRK,h)
```
The solution is returned in an binary format file

### Step 3: Analyzing the solution


#### Energy-Error

First, you must read the solution using **Readbin()** function and convert it to barycentric coordinates:
```julia
tt,U=Readbin(nout,neq,sol1)
solu=ChangeHeltoBar(solU,Gm,tt,nbodyH);
```

```julia
setprecision(BigFloat, 256)

E0=NBodyHam(neq,BigFloat.(u[1]),BigFloat.(Gm))
ΔE = map(x->NBodyHam(neq,BigFloat.(x),BigFloat.(Gm)), u)./E0.-1

plot(tt,log10.(abs.(ΔE)))

```

![Energy plot](/Energy.png)

## More Examples

[Examples](https://github.com/mikelehu/FCIRK/tree/master/Examples)

## References

- [Antoñana, M., Makazaga, J., Murua, Ander. **New Integration Methods for Perturbed ODEs Based
on Symplectic Implicit Runge–Kutta Schemes
with Application to Solar System Simulation**  J Sci Comput 2018](https://doi.org/10.1007/s10915-017-0634-1)

- [Antoñana, M., Makazaga, J., Murua, Ander. **Implicit symplectic methods for high-precision long term integrations of the Sola System**  ??? 2021](https://doi.org/10.1007/s10915-017-0634-1)
