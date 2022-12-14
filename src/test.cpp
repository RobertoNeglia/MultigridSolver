#include "Matrix/dense_matrix.hpp"
#include "Matrix/sparse_matrix_r.hpp"
#include "Poisson2D.hpp"
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
generate_discretized_matrix(SparseMatrix &A, int m, int n) {
  int mn = m * n;
    for (int i = 0; i < mn; i++) {
      A.insert_coeff(4, i, i);
      if (i > 0 && i % m != 0)
        A.insert_coeff(-1, i, i - 1);
      if (i > m - 1)
        A.insert_coeff(-1, i, i - m);
      if (i < mn - 1 && (i % m) != (m - 1))
        A.insert_coeff(-1, i, i + 1);
      if (i < mn - m)
        A.insert_coeff(-1, i, i + m);
    }
}

int
main() {
  Poisson2D poisson(0.0, 1.5, 0.0, 2.0, 1.0 / 2.0);
  poisson.setup();
  poisson.solve();

  // MatrixXi pippo(3, 3);
  //   for (auto i = 0; i < pippo.rows(); i++) {
  //       for (auto j = 0; j < pippo.cols(); j++) {
  //         pippo.insert_coeff(i + j, i, j);
  //       }
  //   }

  //   for (auto i = 0; i < pippo.rows(); i++) {
  //       for (auto j = 0; j < pippo.cols(); j++) {
  //         std::cout << pippo.coeff(i, j) << "\t";
  //       }
  //     std::cout << std::endl;
  //   }

  // unsigned int k = 0;
  // k--;
  // std::cout << k << std::endl;

  //   Domain2D D;
  //   D.initialize(0.0, 1.0, 0.0, 1.5, 1.0 / 2.0);
  //   D.print_domain();
  //   unsigned int mn = D.rows() * D.cols();

  // std::cout << "rows: " << D.rows() << " - cols: " << D.cols() << std::endl;
  // std::cout << "Domain points: " << std::endl;
  // D.print_domain(std::cout);

  // std::cout << "DOMAIN: OK" << std::endl;

  // unsigned int mn = D.rows() * D.cols();
  // std::cout << "rows * cols: " << mn << std::endl;

  //   SparseMatrix A;
  //   A.initialize(mn, mn);
  //   generate_discretized_matrix(A, D.rows(), D.cols());
  //   A.print_matrix();

  // AlgebraicMultigrid amg(D, A);
  return 0;
}