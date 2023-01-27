# AlgebraicMultiGrid

## HOW TO COMPILE AND BUILD:

starting from the root of the repository:
- mkdir build 
- cd build 
- cmake ..
- make

## HOW TO EXECUTE

after compiling, 4 executables will be created inside the build folder:

- SERIAL_JACOBI for a serial implementation of Jacobi
- PARALLEL_JACOBI for a parallel implementation of Jacobi
- SERIAL_MG_METHODS for a serial implementation of both the multigrid methods (geometric and algebraic version)
- PARALLEL_MG_METHODS for a parallel implementation of both the multigrid methods

all the executables can accept a parameter when executed, which corresponds to the size of the coefficient matrix of the linear system

no parameter will launch the demo with a matrix size of 10

### EXAMPLE: ./SERIAL_JACOBI will solve a linear system with system matrix of size 10x10 with a serial implementation of Jacobi

### EXAMPLE: ./PARALLEL_MG_METHODS 32 will solve a linear system in parallel with system matrix of size 32x32
