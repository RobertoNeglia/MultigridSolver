#include "algebraic_multigrid.hpp"
#include "matrix_r.hpp"

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

int
main() {
  int m = 3;
  int n = 3;

  Domain omega;

  Matrix M(m, n);
  M.print_matrix();
  M.print_structure();
  std::cout << std::endl;

  M.insert_coeff(1, 0, 0);
  M.print_matrix();
  M.print_structure();
  std::cout << std::endl;

  M.insert_coeff(1, 1, 1);
  M.print_matrix();
  M.print_structure();
  std::cout << std::endl;

  // SparseMatrix P(m, n);
  // int          mn = m * n;
  // int          k  = 0;
  // for (int i = 0; i < m; i++) {
  //               for (int j = 0; j < n; j++) {
  //                 P.coeffRef(i, j) = k;
  //                 k++;
  //                 std::cout << P.coeff(i, j) << " ";
  //     }
  //   std::cout << std::endl;
  // }
  // SparseMatrix A(mn, mn);

  //   for (int i = 0; i < mn; i++) {
  //     A.coeffRef(i, i) = 4;
  //     if (i > 0 && i % m != 0)
  //       A.coeffRef(i, i - 1) = -1;
  //     if (i > m - 1)
  //       A.coeffRef(i, i - m) = -1;
  //     if (i < mn - 1 && (i % m) != (m - 1))
  //       A.coeffRef(i, i + 1) = -1;
  //     if (i < mn - m)
  //       A.coeffRef(i, i + m) = -1;
  //   }

  // std::cout << A << std::endl;

  // AlgebraicMultigrid amg(A);
  return 0;
}