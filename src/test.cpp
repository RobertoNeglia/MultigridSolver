#include "Matrix/dense_matrix.hpp"
#include "Matrix/dynamic_dense_matrix.hpp"
#include "Matrix/sparse_matrix_r.hpp"
#include "Poisson2D.hpp"
#include "algebraic_multigrid.hpp"
#include "chrono"
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

void
print_vector(std::vector<double> v) {
    for (auto i : v) {
      std::cout << i << " - ";
    }
  std::cout << std::endl;
}

bool
equal_to(std::vector<double> v, double val) {
  bool eq = true;
    for (unsigned int i = 0; i < v.size(); i++) {
      if (std::abs(val - v[i]) > 1.e-6)
        eq = false;
    }

  return eq;
}

auto
timeit(const std::function<void()> &f) {
  using namespace std::chrono;
  const auto start = high_resolution_clock::now();
  f();
  const auto end = high_resolution_clock::now();
  return duration_cast<milliseconds>(end - start).count();
}

int
main() {
  // unsigned int dim = 2;
  // SparseMatrix N;
  // N.initialize(dim * dim, dim * dim);

  // N.insert_coeff(10, 0, 3);
  // N.insert_coeff(10, 1, 2);
  // N.insert_coeff(10, 2, 3);
  // N.insert_coeff(10, 0, 1);
  // N.insert_coeff(10, 1, 0);
  // N.insert_coeff(10, 3, 1);

  // N.print_matrix();

  // SparseMatrix M = N.transpose();

  // M.print_matrix();

  // 4x3
  // Poisson2D poisson(0.0, 1.5, 0.0, 2.0, 1.0 / 2.0);

  // 4x5
  // Poisson2D poisson(0.0, 1.25, 0.0, 1.0, 1 / 4.0);

  // 5x5
  // Poisson2D poisson(0.0, 1.0, 0.0, 1.0, 1 / 5.0);

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

  // exact sol
  std::vector<double> x(A.cols(), 10.);
  // system rhs
  std::vector<double> b = A.mul(x).first;
  // Jacobi initial guess
  std::vector<double> x_guess_jac(A.cols(), 1.);

  unsigned int jac_iter = 5000;

  Jacobi jac(A, b, 1.e-10, jac_iter);

  auto dt = timeit([&]() { jac.solve(x_guess_jac); });
  // jac.solve(x_guess_jac);
  std::cout << "NAIVE JAC: time elapsed " << dt << " [ms]" << std::endl;
  std::cout << "NAIVE JAC: n_iter " << jac.get_iter() << std::endl;
  std::cout << "NAIVE JAC: tol_achived " << jac.get_tol_achieved() << std::endl;

  std::vector<double> e_jac = Jacobi::subvec(x, x_guess_jac);
  // print_vector(e_jac);

  // AMG initial guess
  std::vector<double> x_guess(A.cols(), 1.);

  AlgebraicMultigrid amg(A, b, 20, 20, 1.e-10, 500);

  int flag;
  dt = timeit([&]() { flag = amg.solve(x_guess); });
  // flag = amg.solve(x_guess /*,x*/);

  std::cout << "AMG TIME ELAPSED: " << dt << " [ms]" << std::endl;
  std::cout << "AMG FLAG: " << flag << std::endl;
  std::cout << "AMG TOT ITERATIONS: " << amg.get_iter() << std::endl;
  std::cout << "AMG TOLERANCE ACHIEVED: " << amg.get_tol_achieved() << std::endl;
  std::cout << "AMG JACOBI TOT_ITER: " << amg.get_jac_tot_iter() << std::endl;

  if (equal_to(x_guess, 10.0))
    std::cout << "AMG correct solution found" << std::endl;
  else
    std::cout << ":(" << std::endl;

  // if (equal_to(x_guess_jac, 10.0))
  //   std::cout << "JAC correct solution found" << std::endl;
  // else
  //   std::cout << ":(" << std::endl;

  return 0;
}