#include <iostream>

#include "algebraic_multigrid.hpp"
#include <Eigen/Eigen>

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

using Matrix = Eigen::SparseMatrix<double>;

int
main() {
  int m = 4;
  int n = 4;

  Matrix P(m, n);
  int    mn = m * n;
  int    k  = 0;
    for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
              P.coeffRef(i, j) = k;
              k++;
              std::cout << P.coeff(i, j) << " ";
        }
      std::cout << std::endl;
    }
  Matrix A(mn, mn);

    for (int i = 0; i < mn; i++) {
      A.coeffRef(i, i) = 4;
      if (i > 0 && i % m != 0)
        A.coeffRef(i, i - 1) = -1;
      if (i > m - 1)
        A.coeffRef(i, i - m) = -1;
      if (i < mn - 1 && (i % m) != (m - 1))
        A.coeffRef(i, i + 1) = -1;
      if (i < mn - m)
        A.coeffRef(i, i + m) = -1;
    }

  std::cout << A << std::endl;

  AlgebraicMultigrid amg(A);
  return 0;
}