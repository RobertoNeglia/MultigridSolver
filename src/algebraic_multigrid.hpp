#include <float.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <set>

#include <Eigen/Eigen>

using Matrix = Eigen::SparseMatrix<double>;

// -5 -10 5 2 45 -23 0
// 5 10 -5 -2 -45 23

class AlgebraicMultigrid {
private:
  Matrix                       A;
  std::map<int, std::set<int>> strength;
  const double                 alpha = 1 / 4;

  void
  construct_strength_set() {
      for (int i = 0; i < A.rows(); i++) {
        double max = -DBL_MAX + 10;
          for (int k = 0; k < A.cols(); k++) {
            double a_ik = A.coeff(i, k);
            if (-a_ik > max)
              max = -a_ik;
          }

          for (int j = 0; j < A.cols(); j++) {
            double a_ij = A.coeff(i, j);
            if (j != i && -a_ij >= alpha * max)
              strength.at(i).insert(j);
          }
      }
  }

public:
  AlgebraicMultigrid(Matrix A) : A(A) {}
};