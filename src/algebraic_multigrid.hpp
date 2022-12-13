#ifndef __ALGEBRAIC_MULTIGRID_H_
#define __ALGEBRAIC_MULTIGRID_H_

#define DEBUG 0

#include <float.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <set>

#include "domain.hpp"
#include <Eigen/Eigen>

class AlgebraicMultigrid {
private:
  // Alpha parameter that determines strong connection
  const double alpha = 1.0 / 4.0;

  // Domain of the problem
  Domain2D D;

  // Matrix A of the linear system
  Eigen::SparseMatrix<double> A;

  // Strength set - to each point i is associated the set of points j that influence i
  //                  SET OF POINTS j THAT INFLUENCE i
  std::map<int, std::set<unsigned int>> map_point_to_strength_set;

  // Strength transpose set - to each point i is associated the set of points j that are
  // influenced by i
  //                  SET OF POINTS j THAT ARE INFLUENCED BY i
  std::map<int, std::set<unsigned int>> map_point_to_strength_transpose_set;

  // Associates to each point i the number of points that are influenced by i
  std::map<unsigned int, unsigned int> map_point_to_number_of_influences;

  unsigned int max_lambda;

  Eigen::MatrixXi strength_matrix;

  std::set<unsigned int> C_points;
  std::set<unsigned int> F_points;

  std::pair<unsigned int, unsigned int>
  max_min(unsigned int a, unsigned int b) {
    return a < b ? std::make_pair(b, a) : std::make_pair(a, b);
  }

  double
  find_negative_max_in_row(const unsigned int row) const {
    double max = -DBL_MAX + 10;
      for (unsigned int k = 0; k < A.cols(); k++) {
        double a_ik = A.coeff(row, k);
        if (k != row && -a_ik > max)
          max = -a_ik;
      }

    return max;
  }

  void
  construct_strength_sets() {
    // Constructs the strength set for each point i
      for (int i = 0; i < A.rows(); i++) {
        // finds the max negative coefficient on the row, except diagonal
        double max = find_negative_max_in_row(i);

          for (int j = 0; j < A.cols(); j++) {
            double a_ij = A.coeff(i, j);
              if (j != i && -a_ij >= alpha * max) {
                // insert j to the set of i - j influences i
                map_point_to_strength_set[i].insert(j);
            }
          }
      }
    // Constructs the strength transpose set for each point i
      for (auto [j, s_j] : map_point_to_strength_set) {
        // for each point j
          for (auto i : s_j) {
            // add j to the strength transpose set of i, each i is inside the strength set
            // of j
            map_point_to_strength_transpose_set[i].insert(j);
          }
      }

      for (auto [i, s_i] : map_point_to_strength_transpose_set) {
        map_point_to_number_of_influences[i] = s_i.size();
      }
  }

  void
  construct_strength_matrix() {
    unsigned int max = 0;
      for (auto [i, lambda_i] : map_point_to_number_of_influences) {
        std::pair<unsigned int, unsigned int> coordinates               = D.onedim2twodim(i);
        strength_matrix.coeffRef(coordinates.first, coordinates.second) = lambda_i;
        if (lambda_i > max)
          max = lambda_i;
      }

    max_lambda = max;
  }

  bool
  inside_boundary(const unsigned int i, const unsigned int j) const {
    // no need of i >= 0 and j >= 0 checks since they're unsgined
    return i < D.rows() && j < D.cols();
  }

  bool
  not_C_point(const unsigned int i, const unsigned int j) const {
    return (C_points.find(D.twodim2onedim(i, j)) == C_points.end());
  }

  bool
  not_F_point(const unsigned int i, const unsigned int j) const {
    return (F_points.find(D.twodim2onedim(i, j)) == F_points.end());
  }

  void
  update_strength(const unsigned int i, const unsigned int j) {
      if (inside_boundary(i, j) && not_C_point(i, j) && not_F_point(i, j)) {
        unsigned int &lambda =
          reinterpret_cast<unsigned int &>(strength_matrix.coeffRef(i, j));
        lambda++;
        if (lambda > max_lambda)
          max_lambda = lambda;
    }
  }

  void
  update_neighbours_strength(const unsigned int i, const unsigned int j) {
    update_strength(i - 1, j);
    update_strength(i, j - 1);
    update_strength(i + 1, j);
    update_strength(i, j + 1);

    //   if (inside_boundary(i - 1, j) && not_C_point(i - 1, j) && not_F_point(i - 1, j)) {
    //     strength_matrix.coeffRef(i - 1, j)++;
    // }
    //   if (inside_boundary(i, j - 1) && not_C_point(i, j - 1) && not_F_point(i, j - 1)) {
    //     strength_matrix.coeffRef(i, j - 1)++;
    // }
    //   if (inside_boundary(i + 1, j) && not_C_point(i + 1, j) && not_F_point(i + 1, j)) {
    //     strength_matrix.coeffRef(i + 1, j)++;
    // }
    //   if (inside_boundary(i, j + 1) && not_C_point(i, j + 1) && not_F_point(i, j + 1)) {
    //     strength_matrix.coeffRef(i, j + 1)++;
    // }
  }

  void
  make_F_point(const unsigned int i, const unsigned int j) {
      if (inside_boundary(i, j) && not_C_point(i, j) && not_F_point(i, j)) {
        F_points.insert(D.twodim2onedim(i, j));
        update_neighbours_strength(i, j);
    }
  }

  void
  neighbours_to_F_points(const unsigned int i, const unsigned int j) {
    make_F_point(i - 1, j);
    make_F_point(i, j - 1);
    make_F_point(i + 1, j);
    make_F_point(i, j + 1);

    //   if (inside_boundary(i - 1, j) && not_C_point(i - 1, j)) {
    //     F_points.insert(D.twodim2onedim(i - 1, j));
    //     std::cout << "UP: OK" << std::endl;
    //     update_neighbours_strength(i - 1, j);
    // }
    //   if (inside_boundary(i, j - 1) && not_C_point(i, j - 1)) {
    //     F_points.insert(D.twodim2onedim(i, j - 1));
    //     std::cout << "LEFT: OK" << std::endl;
    //     update_neighbours_strength(i, j - 1);
    // }
    //   if (inside_boundary(i + 1, j) && not_C_point(i + 1, j)) {
    //     F_points.insert(D.twodim2onedim(i + 1, j));
    //     std::cout << "DOWN: OK" << std::endl;
    //     update_neighbours_strength(i + 1, j);
    // }
    //   if (inside_boundary(i, j + 1) && not_C_point(i, j + 1)) {
    //     F_points.insert(D.twodim2onedim(i, j + 1));
    //     std::cout << "RIGHT: OK" << std::endl;
    //     update_neighbours_strength(i, j + 1);
    // }
  }

  void
  make_C_point(unsigned int i, unsigned int j) {
    C_points.insert(D.twodim2onedim(i, j));
    neighbours_to_F_points(i, j);
      if (DEBUG) {
        print_C_points();
        print_F_points();
        print_strength_matrix();
        press_to_continue();
    }
  }

  void
  coarsen() {
      for (int i = strength_matrix.rows() - 1; i > -1; i--) {
          for (int j = 0; j < strength_matrix.cols(); j++) {
            unsigned int lambda = strength_matrix.coeff(i, j);
              if (lambda >= max_lambda && not_C_point(i, j) && not_F_point(i, j)) {
                make_C_point(i, j);
            }
          }
      }
  }

  void
  print_strength_sets(std::ostream &os = std::cout) const {
      for (auto [i, s_i] : map_point_to_strength_set) {
        os << "Points which influence point " << i << ":" << std::endl;

          for (auto j : s_i) {
            os << "\t" << j << std::endl;
          }
      }

      for (auto [i, s_i] : map_point_to_strength_transpose_set) {
        os << "Points that are influenced by point " << i << ":" << std::endl;

          for (auto j : s_i) {
            os << "\t" << j << std::endl;
          }
      }

      for (int i = 0; i < A.rows(); i++) {
        os << map_point_to_number_of_influences.at(i);
      }
    os << std::endl;

      for (auto [i, lambda_i] : map_point_to_number_of_influences) {
        os << "Number of points influenced by " << i << ": " << lambda_i << std::endl;
      }
  }

  void
  print_strength_matrix(std::ostream &os = std::cout) const {
    os << "Strenght matrix is: " << std::endl;
      for (int i = 0; i < strength_matrix.rows(); i++) {
          for (int j = 0; j < strength_matrix.cols(); j++) {
            os << strength_matrix.coeff(i, j) << "\t";
          }
        os << std::endl;
      }
  }

  void
  print_C_points(std::ostream &os = std::cout) const {
    os << "C-points are: " << std::endl;
      for (auto i : C_points) {
        os << i << " - ";
      }
    os << std::endl;
  }

  void
  print_F_points(std::ostream &os = std::cout) const {
    os << "F-points are: " << std::endl;
      for (auto i : F_points) {
        os << i << " - ";
      }
    os << std::endl;
  }

  void
  press_to_continue(std::ostream &os = std::cout, std::istream &is = std::cin) const {
    os << "Press enter to continue: ";
    is.ignore();
  }

public:
  AlgebraicMultigrid(Domain2D D, Eigen::SparseMatrix<double> A) : D(D), A(A) {
    strength_matrix.resize(D.rows(), D.cols());
    construct_strength_sets();
    construct_strength_matrix();
      if (DEBUG) {
        print_strength_sets();
        print_strength_matrix();
    }

    coarsen();

    print_strength_matrix();
    print_C_points();
    print_F_points();
  }
};

#endif