#ifndef __ALGEBRAIC_MULTIGRID_H_
#define __ALGEBRAIC_MULTIGRID_H_

#include <float.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <set>

#include <Eigen/Eigen>

// -5 -10 5 2 45 -23 0
// 5 10 -5 -2 -45 23

class AlgebraicMultigrid {
private:
  Eigen::SparseMatrix<double>  A;
  std::map<int, std::set<int>> strength;
  std::map<int, std::set<int>> strength_transpose;
  std::map<int, int>           number_of_influences;
  const double                 alpha = 1.0 / 4.0;

  void
  construct_strength_sets() {
      for (int i = 0; i < A.rows(); i++) {
        double max = -DBL_MAX + 10;

          for (int k = 0; k < A.cols(); k++) {
            double a_ik = A.coeff(i, k);
            if (k != i && -a_ik > max)
              max = -a_ik;
          }

        // std::cout << "Max of row " << i << " is: " << max << std::endl;

          for (int j = 0; j < A.cols(); j++) {
            double a_ij = A.coeff(i, j);
              if (j != i && -a_ij >= alpha * max) {
                // std::cout << "a_ij: " << a_ij << std::endl;
                strength[i].insert(j);
            }
          }
      }

      for (auto [j, s_j] : strength) {
          for (auto i : s_j) {
            strength_transpose[i].insert(j);
          }
      }

      for (auto [i, s_i] : strength_transpose) {
        number_of_influences[i] = s_i.size();
      }
  }

  void
  print_strength_sets() {
      for (auto [i, s_i] : strength) {
        std::cout << "Points which influence point " << i << ":" << std::endl;

          for (auto j : s_i) {
            std::cout << "\t" << j << std::endl;
          }
      }

      for (auto [i, s_i] : strength_transpose) {
        std::cout << "Points that are influenced by point " << i << ":" << std::endl;

          for (auto j : s_i) {
            std::cout << "\t" << j << std::endl;
          }
      }

      for (int i = 0; i < A.rows(); i++) {
        std::cout << number_of_influences.at(i);
      }
    std::cout << std::endl;

      for (auto [i, lambda_i] : number_of_influences) {
        std::cout << "Number of points influenced by " << i << ": " << lambda_i
                  << std::endl;
      }
  }

public:
  AlgebraicMultigrid(Eigen::SparseMatrix<double> A) : A(A) {
    construct_strength_sets();
    print_strength_sets();
  }
};

#endif