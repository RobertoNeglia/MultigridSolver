#include "Matrix/dense_matrix.hpp"
#include "Matrix/sparse_matrix_r.hpp"
#include "Poisson2D.hpp"
#include "algebraic_multigrid.hpp"
#include "domain.hpp"
#include "jacobi_r.hpp"

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
  // 4x3
  // Poisson2D poisson(0.0, 1.5, 0.0, 2.0, 1.0 / 2.0);

  // 4x5
  // Poisson2D poisson(0.0, 1.25, 0.0, 1.0, 1 / 4.0);

  // 5x5
  Poisson2D poisson(0.0, 1.0, 0.0, 1.0, 1 / 5.0);

  // 4x4
  // Poisson2D poisson(0.0, 1.0, 0.0, 1.0, 1 / 4.0);

  // 3x3
  // Poisson2D poisson(0.0, 1.5, 0.0, 1.5, 1.0 / 2.0);
  // poisson.setup();
  // poisson.solve();

  Domain2D D;
  D.initialize(0.0, 1.0, 0.0, 1.0, 1.0 / 16.0);

  int          n = D.cols() * D.rows();
  SparseMatrix A;
  A.initialize(n, n);

  generate_discretized_matrix(A, D.rows(), D.cols());
  // for (int i = 0; i < n; i++) {
  //   A.insert_coeff(2, i, i);
  //   if (i > 0)
  //     A.insert_coeff(-1, i, i - 1);
  //   if (i < n - 1)
  //     A.insert_coeff(-1, i, i + 1);
  // }

  // exact sol
  std::vector<double> x(A.cols(), 10.);
  // initial guess
  std::vector<double> x_guess(A.cols(), -100.);

  std::vector<double> b = A.mul(x).first;

  AlgebraicMultigrid amg(D, A, b);

  amg.solve(x_guess);

  Jacobi jac(A, b, 1.e-6, 1000);

  std::vector<double> x_guess_jac(A.cols(), -100.);

  jac.solve(x_guess_jac);

  std::cout << "jac tot iter: " << jac.get_iter();

  // A.print_matrix();

  // unsigned int max_it = 1000;
  // double       tol    = 1.e-9;

  // Jacobi solver(A, b, tol, max_it);

  // solver.solve(x_guess);

  // std::cout << "n_it - " << max_it << std::endl;
  // std::cout << "tol - " << tol << std::endl;

  //   for (auto i : x_guess) {
  //     std::cout << i << std::endl;
  //   }

  // int    n = 10;
  // Matrix A(n, n);

  //   for (int i = 0; i < n; i++) {
  //     A.insertElement(2, i, i);
  //     if (i > 0)
  //       A.insertElement(-1, i, i - 1);
  //     if (i < n - 1)
  //       A.insertElement(-1, i, i + 1);
  //   }

  // A.print();

  // std::vector<double> x(A.getCols(), 10.);

  // std::vector<double> b = A.Multiplication(x);

  // std::vector<double> x_guess(A.getCols(), 0.);

  // unsigned int max_it = 1000;
  // double       tol    = 1.e-9;

  // Jacobi solver(A, b, tol, max_it);

  // solver.solve(x_guess);

  // std::cout << "n_it" << max_it << std::endl;
  // std::cout << "tol" << tol << std::endl;

  //   for (auto i : x_guess) {
  //     std::cout << i << std::endl;
  //   }

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

  // std::cout << "rows: " << D.rows() << " - cols: " << D.cols() <<
  // std::endl; std::cout << "Domain points: " << std::endl;
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