#include "Matrix/matrix_r.hpp"
#include "algebraic_multigrid.hpp"

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

void
generate_discretized_matrix(Eigen::SparseMatrix<double> &A, int m, int n) {
  int mn = m * n;
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
}

int
main() {
  int m = 3;
  int n = 3;

  Matrix M(m, n);
  M.print_matrix();
  M.print_structure();
  std::cout << std::endl;

  M.insert_coeff(5, 0, 1);
  M.insert_coeff(6, 1, 0);
  M.insert_coeff(8, 1, 1);
  M.insert_coeff(7, 2, 2);
  M.insert_coeff(-30, 0, 2);
  M.insert_coeff(15, 2, 0);
  M.print_matrix();
  M.print_structure();
  std::cout << std::endl;

  std::vector<double> v{9, 3, -1};

  std::cout << "vector v: " << std::endl;
  for (auto i : v)
    std::cout << i << " - ";
  std::cout << std::endl;

  std::pair<std::vector<double>, bool> Mv = M.mul(v);
    if (Mv.second) {
      std::cout << "vector M*v: " << std::endl;
      for (auto i : Mv.first)
        std::cout << i << " - ";
      std::cout << std::endl;
  } else
    std::cout << "ERROR: couldn't compute multiplication" << std::endl;

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
  //   int                         mn = m * n;
  //   Eigen::SparseMatrix<double> A(mn, mn);
  //   generate_discretized_matrix(A, m, n);

  //   std::cout << A << std::endl;

  //   AlgebraicMultigrid amg(A);
  return 0;
}