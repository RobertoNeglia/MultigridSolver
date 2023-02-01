//#include <mpi.h>

#include <iostream>

#include "matrixz.hpp"
int
main() {
  int    n = 3;
  Matrix m = Matrix(n, n);
  m.print();
  m.insertElement(6, 0, 0);
  m.insertElement(7, 2, 2);
  m.insertElement(5, 1, 1);

  m.print();

  m.insertElement(4, 1, 2);

  m.print();
}
/*
    TODO:
        - discretizzazione problema
        - SERIAL
            - implementazione matrice
            - implementazione smoother (Jacobi/GS)
            - coarsening

        - PARALLEL
            - parallelizzare smoother
            - parallelizzare coarsening
*/
