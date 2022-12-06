//#include <mpi.h>


#include <iostream>
#include "matrix.hpp"
int
main() {
    Matrix m = Matrix(3,3);
    m.print();
    m.insertElement(3,0,1);
    m.print();
    m.insertElement(-3,1,1);
    m.insertElement(0,1,1);
    m.insertElement(1,0,2);
    //m.insert
    m.insertElement(1,2,0);
    m.insertElement(-1,2,2);
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

