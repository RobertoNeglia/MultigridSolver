#include <Eigen/Eigen>

int
main() {
  Eigen::SparseMatrix<double> M(10, 10);

  M.transpose();
}