# Geometric Multigrid solver project

## Partecipants

Beltrame Leonardo \
Davide Zaira \
Neglia Roberto 

## Short description

Given the number of elements N (**assumed to be even**) to split the domain $\Omega = (0, 1)^2$ into, the system matrix $A$ of size $n^2\times n^2$ is assembled, $n = N - 1$,  with the following structure that arise from the discretization of the 2D Poisson problem:

![system_matrix](/imgs/system_matrix.png)

with $D$, of size $n\times n$, defined as 

![D_matrix](/imgs/D_matrix.png)

and $I$, also of size $n\times n$, is the identity matrix.

The iterative solver solves the linear system $Ax = b$.

The matrix is stored in a CSR (Compressed Sparse Row) format to exploit the sparsity that arise from the problem discretization.

The code is completely general and can be used with any symmetric positive definite matrix (to assure convergence, from iterative solvers theory).

The smoother used is Jacobi, but the code is written so that it can be changed to Gauss-Seidel or any other iterative solver that inherits from the `IterativeSolver` abstract class.

the `GeometricMultigrid` class is the implementation of the 2 level multigrid solver, while the `MultilevelGeometricMultigrid` is the implementation of the multilevel multigrid solver, which should be used for a number of levels **greater or equal than 3**: using this class for a 2 level multigrid solver will result in a **non-efficient solver** due to the particular data structures used for the multilevel implementation.

There is also an implementation of the `AlgebraicMultigrid` solver: it works, does the CF-splitting correctly, but it's inefficient, probably due to the kind of interpolation we implemented.

OpenMP has been used to parallelize the code, in particular all matrix-vector multiplications and vector algebra operations.


## How to compule and build:

Starting from the root of the repository:
```
mkdir build 
cd build 
cmake ..
make
```

## How to execute

After compiling, 6 executables will be created inside the build folder:

- `SERIAL_JACOBI` for a serial implementation of Jacobi;
- `PARALLEL_JACOBI` for a parallel implementation of Jacobi;
- `SERIAL_MG_METHODS` for a serial implementation of both the two level multigrid methods (geometric and algebraic version);
- `PARALLEL_MG_METHODS` for a parallel implementation of both the two level multigrid methods;
- `SERIAL_ML_MG_METHODS` for a serial implementation of the multilevel multigrid method;
- `PARALLEL_ML_MG_METHODS` for a parallel implementation of the multilevel multigrid method.

All the executables can accept a parameter when executed, which corresponds to the N above (system matrix will have size $n^2\times n^2$).

No parameter will launch the demo with a matrix size of $N = 10$ (matrix size $81\times 81$).

For the multilevel multigrid demo, a second parameter can be passed to the program which corresponds to the number of levels of multigrid: no second parameter will launch the multilevel multigrid demo with 3 levels.

## Examples

- EXAMPLE 1: `./SERIAL_JACOBI` will solve a linear system with system matrix of size $81\times 81$ with a serial implementation of Jacobi

- EXAMPLE: `./PARALLEL_MG_METHODS 32` will solve a linear system in parallel with system matrix of size $961\times 961$

- EXAMPLE: `./PARALLEL_ML_MG_METHODS 64 4` will solve a linear system in parallel with system matrix of size $3969\times 3969$
