#include <vector>

#include "matrix.hpp"

class Jacobi {
public:
  Jacobi(Matrix A, std::vector<double> b) : A(A), b(b) {
    tol = 1.e-6;
  }

  Jacobi(Matrix A, std::vector<double> b, double tol) : A(A), b(b), tol(tol) {}

  void
  solve(std::vector<double> x0){

  };

  void
  setTol(double tol) {
    tol = tol;
  }

private:
  Matrix              A;
  std::vector<double> b;
  double              tol;
};