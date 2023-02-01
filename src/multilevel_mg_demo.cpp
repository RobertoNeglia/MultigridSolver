#include <chrono>
#include <functional>

#include "LinearAlgebra/IterativeSolvers/multilevel_geometric_multigrid.hpp"

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

auto
timeit(const std::function<void()> &f) {
  using namespace std::chrono;
  const auto start = high_resolution_clock::now();
  f();
  const auto end = high_resolution_clock::now();
  return duration_cast<milliseconds>(end - start).count();
}

int
main(int argc, char **argv) {
  unsigned int N        = 10;
  unsigned int n_levels = 3;
    if (argc > 1) {
      N = std::atoi(argv[1]);
  }
    if (argc > 2) {
      n_levels = std::atoi(argv[2]);
  }

  unsigned int n = (N - 1) * (N - 1);

  SparseMatrix A;
  A.initialize(n, n);

  generate_discretized_matrix(A, N - 1, N - 1);

  double exact_sol = 1.0;
  // exact solution of the linear system
  std::vector<double> x(A.cols(), exact_sol);
  // system rhs
  std::vector<double> b = A.mul(x);
  // print_vector(b);

  // GMG initial guess
  double         initial_guess = 0.0;
  Vector<double> gmg_guess(A.cols(), initial_guess);
  Vector<double> amg_guess(A.cols(), initial_guess);

  const double       tol      = 1.e-8;
  const unsigned int max_iter = 10000;

  const unsigned int mg_pre_nu  = 50;
  const unsigned int mg_post_nu = 50;

  MultilevelGeometricMultigrid gmg(A, b, mg_pre_nu, mg_post_nu, tol, max_iter, n_levels);
  gmg.setup();

  int  flag;
  auto dt = timeit([&]() { flag = gmg.solve(gmg_guess, n_levels); });

  std::cout << "GMG SOLVE TIME ELAPSED: " << dt << " [ms]" << std::endl;
  std::cout << "GMG FLAG: " << flag << std::endl;
  std::cout << "GMG TOT ITERATIONS: " << gmg.get_iter() << std::endl;
  std::cout << "GMG TOLERANCE ACHIEVED: " << gmg.get_tol_achieved() << std::endl;
  std::cout << "GMG JACOBI TOT_ITER: " << gmg.get_tot_smoother_iter() << std::endl;

  if (equal_to(gmg_guess, exact_sol))
    std::cout << "GMG correct solution found" << std::endl;
  else
    std::cout << ":(" << std::endl;

  std::cout << "===============================================" << std::endl << std::endl;

#ifdef _OPENMP
  std::cout << "ER SAMBUCONE MOLINARI CURATIVO" << std::endl;
#endif

  return 0;
}