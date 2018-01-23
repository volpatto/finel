# FINEL -- A general purpose FINite ELement code

I present here my implementation of the Finite Element Method with Galerkin Approximation based on the classical book Finite Elements: An Introduction wrote by Tinsley Oden, Eric Becker and Graham Carey (such a great book!!!). Also, I inspired myself by some ideas showed in The Finite Element Method: Linear Static and Dynamic Finite Element Analysis, Thomas Hughes. The code provided there (DLEARN) was a great help for me.

## Features

* Modular Fortran programming with modern features (f95 and f03);
* Derived Data Structure;
* Reader for Triangle and EasyMesh mesh generators;
* Element types available:
    * Linear, quadratic and cubic to 1D problems;
    * Linear triangles to 2D problems.
* In-house numerical linear algebra solvers available:
    * Jacobi;
    * Gauss-Seidel;
    * SOR;
* Direct linear algebra solvers available:
    * Gauss foward substitution with triangular superior reduction.
* Picard method to solve nonlinear problems;
* Writer result data to the following file extensions:
    * .dat;
    * .csv;
    * .vtk.
* General 1D/2D stationary or transient advection-diffusion-reaction equation;
* Stabilization methods:
    * SUPG;
    * GGLS;
    * CAU (testing stage).

## Contact

My name is Diego, but you can call me Volps. I'm a Doctoral Student at National Laboratory for Scientific Computing (LNCC/MCTIC, Brazil), sited in Petrópolis/RJ. Feel free to send me an e-mail: volpatto@lncc.br.

Greetings!
