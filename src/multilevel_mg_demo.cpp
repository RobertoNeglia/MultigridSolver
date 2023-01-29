#include <chrono>
#include <functional>

#include "LinearAlgebra/IterativeSolvers/algebraic_multigrid.hpp"
#include "LinearAlgebra/IterativeSolvers/geometric_multigrid.hpp"
#include "LinearAlgebra/sparse_matrix.hpp"

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
      if (std::abs(val - v[i]) > 1.e-4)
        eq = false;
    }

  return eq;
}

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
  unsigned int N = 10;
    if (argc == 2) {
      N = std::atoi(argv[1]);
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

  const unsigned int mg_pre_nu  = 150;
  const unsigned int mg_post_nu = 150;

  GeometricMultigrid gmg(A, b, mg_pre_nu, mg_post_nu, tol, max_iter);
  gmg.setup();

  int  flag;
  auto dt = timeit([&]() { flag = gmg.solve(gmg_guess, 3); });

  std::cout << "GMG TIME ELAPSED: " << dt << " [ms]" << std::endl;
  std::cout << "GMG FLAG: " << flag << std::endl;
  std::cout << "GMG TOT ITERATIONS: " << gmg.get_iter() << std::endl;
  std::cout << "GMG TOLERANCE ACHIEVED: " << gmg.get_tol_achieved() << std::endl;
  std::cout << "GMG JACOBI TOT_ITER: " << gmg.get_tot_smoother_iter() << std::endl;

  if (equal_to(gmg_guess, exact_sol))
    std::cout << "GMG correct solution found" << std::endl;
  else
    std::cout << ":(" << std::endl;

  std::cout << "===============================================" << std::endl << std::endl;

  AlgebraicMultigrid amg(A, b, mg_pre_nu, mg_post_nu, tol, max_iter);
  amg.setup();

  dt = timeit([&]() { flag = amg.solve(amg_guess); });

  std::cout << "AMG TIME ELAPSED: " << dt << " [ms]" << std::endl;
  std::cout << "AMG FLAG: " << flag << std::endl;
  std::cout << "AMG TOT ITERATIONS: " << amg.get_iter() << std::endl;
  std::cout << "AMG TOLERANCE ACHIEVED: " << amg.get_tol_achieved() << std::endl;
  std::cout << "AMG JACOBI TOT_ITER: " << amg.get_tot_smoother_iter() << std::endl;

  if (equal_to(amg_guess, exact_sol))
    std::cout << "AMG correct solution found" << std::endl;
  else
    std::cout << ":(" << std::endl;

#ifdef _OPENMP
  std::cout << "ER SAMBUCONE MOLINARI CURATIVO" << std::endl;
#endif

  return 0;
}