#include "Matrix/sparse_matrix_r.hpp"
#include "algebraic_multigrid.hpp"
#include "domain.hpp"

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
  int m = 5;
  int n = 4;

  // Matrix M(m, n);
  // M.print_matrix();
  // M.print_structure();
  // std::cout << std::endl;

  // M.insert_coeff(5, 0, 1);
  // M.insert_coeff(6, 1, 0);
  // M.insert_coeff(8, 1, 1);
  // M.insert_coeff(7, 2, 2);
  // M.insert_coeff(-30, 0, 2);
  // M.insert_coeff(15, 2, 0);
  // M.print_matrix();
  // M.print_structure();
  // std::cout << std::endl;

  // std::vector<double> v{0, 0, 0};

  // std::cout << "vector v: " << std::endl;
  // for (auto i : v)
  //   std::cout << i << " - ";
  // std::cout << std::endl;

  // std::pair<std::vector<double>, bool> Mv = M.mul(v);
  //   if (Mv.second) {
  //     std::cout << "vector M*v: " << std::endl;
  //     for (auto i : Mv.first)
  //       std::cout << i << " - ";
  //     std::cout << std::endl;
  // } else
  //   std::cout << "ERROR: couldn't compute multiplication" << std::endl;

  Domain2D D(m, n);

  std::cout << "Domain points: " << std::endl;
  D.print_domain(std::cout);

  int mn = m * n;

  // unsigned int p = 33;
  // unsigned int i, j;

  // i = p / D.cols();
  // j = p % D.cols();

  // std::cout << "p: " << p << " - (i,j): (" << i << "," << j << ")" << std::endl;

  Eigen::SparseMatrix<double> A(mn, mn);
  generate_discretized_matrix(A, m, n);

  // std::cout << A << std::endl;

  AlgebraicMultigrid amg(D, A);
  return 0;
}