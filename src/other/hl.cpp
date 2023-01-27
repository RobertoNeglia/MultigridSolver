//#include <mpi.h>


#include <iostream>
#include "matrix.hpp"
int
main() {
    Matrix m = Matrix(4,4);
    for (unsigned int k = 0; k < 4; k++)
    {
        m.insertElement(-1, k, k);
    }
    m.print();

    return 0;
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