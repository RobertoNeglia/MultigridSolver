// #include <mpi.h>

#include <chrono>
#include <iostream>
#include <vector>
// #include "matrix.hpp"
#include "jacobi.hpp"

using namespace std;

int
main() {
  Matrix         m = Matrix(50, 50);
  Matrix         d = Matrix(50, 50);
  vector<double> b(50, -1.0);
  vector<double> x(50, 1.0);
    // x[30] = 1.0;
    for (unsigned int k = 0; k < x.size(); k++) {
      m.insertElement(3, k, k);
        if (k < x.size() - 1) {
          m.insertElement(-1, k + 1, k);
          m.insertElement(-1, k, k + 1);
      }
    }
  m.print();
  d = m.JacPrec();
  d.print();

  Jacobi jac = Jacobi(m, b);
  jac.solve(x);
  cout << endl;
    for (unsigned int i = 0; i < x.size(); i++) {
      cout << x[i] << " ";
    }
  cout << endl;

  // vector<double> prova = m.Multiplication(x);
  // for (unsigned int i = 0; i < x.size() ; i++)
  // {
  //     cout << prova[i] << " ";
  // }
  // cout << endl;
  // b = m.Multiplication(x);
  // cout << "Il vettore risultante e' " << endl;
  // for (unsigned int i = 0; i < b.size(); i++)
  // {
  //     cout << b[i] << " ";
  // }
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
