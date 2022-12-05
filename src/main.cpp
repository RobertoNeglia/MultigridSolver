//#include <mpi.h>


#include <iostream>
#include "matrix.hpp"
int
main() {
    Matrix m = Matrix(3,3);
    m.print();
    m.insertElement(3,0,1);
    m.insertElement(-3,1,1);
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

