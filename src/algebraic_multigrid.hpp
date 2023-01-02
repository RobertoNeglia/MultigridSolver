#ifndef __ALGEBRAIC_MULTIGRID_H_
#define __ALGEBRAIC_MULTIGRID_H_

#define DEBUG 0
#define DBL_MIN -100000000

#include <algorithm>
#include <iostream>
#include <map>
#include <set>

#include "Matrix/dense_matrix.hpp"
#include "Matrix/dynamic_dense_matrix.hpp"
#include "Matrix/sparse_matrix_r.hpp"
#include "domain.hpp"
#include "jacobi_r.hpp"

class AlgebraicMultigrid {
  //---------------------------------------------------------------------------------
  // PRIVATE MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
private:
  /**
   * Alpha parameter to determine strong connection
   */
  const double alpha = 1.0 / 4.0;

  /**
   * Pre&Post smoother number of iterations
   */
  const double nu = 200;

  /**
   * System matrix
   */
  SparseMatrix A;

  /**
   * System right hand side
   */
  std::vector<double> b;

  /**
   * Strength set - to each point i is associated the set of points j that influence i
   *        SET OF POINTS j THAT INFLUENCE POINT i
   *        SET OF POINTS j ON WHICH POINT i DEPENDS
   */
  std::map<int, std::set<unsigned int>> map_point_to_strength_set;

  /**
   * Strength transpose set - to each point i is associated the set of points j that are
   * influenced by i
   *        SET OF POINTS j THAT ARE INFLUENCED BY POINT i
   *        SET OF POINTS j THAT DEPEND ON POINT i
   */
  std::map<int, std::set<unsigned int>> map_point_to_strength_transpose_set;

  // Associates to each point i the number of points that are influenced by i
  std::map<unsigned int, unsigned int> map_point_to_number_of_influences;

  unsigned int max_lambda;

  MatrixI  strength_matrix;
  MatrixDI coarsened_domain;

  std::set<unsigned int> C_points;
  std::set<unsigned int> F_points;

  // Associates to each point i the c-points that are used to interpolate i
  std::map<unsigned int, std::set<unsigned int>> map_point_to_interpolating_points;

  double
  find_negative_max_in_row(const unsigned int row) const {
    double max = DBL_MIN;
      for (unsigned int k = 0; k < A.cols(); k++) {
        double a_ik = A.coeff(row, k).first;
        if (k != row && -a_ik > max)
          max = -a_ik;
      }

    return max;
  }

  void
  construct_strength_sets() {
    // Constructs the strength set for each point i
      for (unsigned int i = 0; i < A.rows(); i++) {
        // finds the max negative coefficient on the row, except diagonal
        double max = find_negative_max_in_row(i);

          for (unsigned int j = 0; j < A.cols(); j++) {
            double a_ij = A.coeff(i, j).first;
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
        if (s_i.size() > max_lambda)
          max_lambda = s_i.size();
      }
  }

  bool
  not_C_point(const unsigned int p) const {
    return !C_points.contains(p);
  }

  bool
  not_F_point(const unsigned int p) const {
    return !F_points.contains(p);
  }

  void
  update_strength(const unsigned int p) {
      if (not_C_point(p) && not_F_point(p)) {
        map_point_to_number_of_influences[p]++;
        unsigned int lambda = map_point_to_number_of_influences[p];
        lambda++;
        if (lambda > max_lambda)
          max_lambda = lambda;
    }
  }

  void
  update_neighbours_strength(const unsigned int p) {
      for (auto i : map_point_to_strength_set[p]) {
        update_strength(i);
      }
  }

  void
  make_F_point(const unsigned int p) {
      if (not_C_point(p) && not_F_point(p)) {
        F_points.insert(p);
        update_neighbours_strength(p);
    }
  }

  void
  neighbours_to_F_points(const unsigned int p) {
      for (auto j : map_point_to_strength_transpose_set[p]) {
        make_F_point(j);
      }
  }

  void
  make_C_point(const unsigned int p) {
    C_points.insert(p);
    neighbours_to_F_points(p);
      if (DEBUG) {
        print_C_points();
        print_F_points();
        // press_to_continue();
    }
  }

  void
  C_F_splitting() {
    const unsigned int n_elements = A.rows();
      while (C_points.size() + F_points.size() < n_elements) {
          for (auto [i, lambda_i] : map_point_to_number_of_influences) {
              if (lambda_i >= max_lambda && not_C_point(i) && not_F_point(i)) {
                // C point found
                make_C_point(i);
            }
          }
        if (max_lambda > 0)
          max_lambda--;
      }
  }

  std::vector<double>
  coarsen_residual(const std::vector<double> &r) {
    std::vector<double> rc(C_points.size());

    double val;
    auto   c_it = C_points.begin();
      for (unsigned int i = 0; i < rc.size(); i++) {
          if (i == 0) {
            val = (r[*c_it] + r[*c_it + 1]) / 2;
          } else if (i == rc.size() - 1) {
            val = (r[*c_it - 1] + r[*c_it]) / 2;
          } else {
            val = (r[*c_it - 1] + r[*c_it] + r[*c_it + 1]) / 3;
          }
        rc[i] = val;
        c_it++;
      }

    return rc;
  }

  SparseMatrix
  coarsen_matrix(const SparseMatrix &A) const {
    SparseMatrix Ac;
    Ac.initialize(C_points.size(), C_points.size());

      for (unsigned int i = 0; i < Ac.rows(); i++) {
          for (unsigned int j = 0; j < Ac.cols(); j++) {
            unsigned int ii = 2 * i, jj = 2 * j;
            double       new_coeff = A.coeff(ii, jj).first;
              // Boundary management
              if (ii == 0) { // elements on first row : only use next row value
                new_coeff += 0.25 * A.coeff(ii + 1, jj).first;
              } else if (ii == A.rows() - 1) { // elements on last row: only use prev row value
                new_coeff += 0.25 * A.coeff(ii - 1, jj).first;
              } else { // otherwise use both
                new_coeff += 0.25 * (A.coeff(ii - 1, jj).first + A.coeff(ii + 1, jj).first);
              }

              if (jj == 0) { // elements on first col : only use next row value
                new_coeff += 0.25 * A.coeff(ii, jj + 1).first;
              } else if (jj == A.cols() - 1) { // elements on last col: only use prev row value
                new_coeff += 0.25 * A.coeff(ii, jj - 1).first;
              } else { // otherwise use both
                new_coeff += 0.25 * (A.coeff(ii, jj - 1).first + A.coeff(ii, jj + 1).first);
              }

            Ac.insert_coeff(new_coeff, i, j);
          }
      }
    return Ac;
  }

  std::vector<double>
  interpolate_error(const std::vector<double> e) {
    std::vector<double> ef(C_points.size() + F_points.size());

      if (*(C_points.begin()) == 0) {
          for (unsigned int i = 0; i < e.size() - 1; i++) {
            ef[2 * i]     = e[i];
            ef[2 * i + 1] = (e[i] + e[i + 1]) / 2;
          }
        ef[ef.size() - 1] = e[e.size() - 1];
      } else {
        ef[0] = e[0];
        ef[1] = e[0];
          for (unsigned int i = 1; i < e.size(); i++) {
            ef[2 * i] = (e[i - 1] + e[i]) / 2;
            if (2 * i + 1 < ef.size())
              ef[2 * i + 1] = e[i];
          }
      }

    return ef;
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

      for (unsigned int i = 0; i < A.rows(); i++) {
        os << map_point_to_number_of_influences.at(i);
      }
    os << std::endl;

      for (auto [i, lambda_i] : map_point_to_number_of_influences) {
        os << "Number of points influenced by " << i << ": " << lambda_i << std::endl;
      }
  }

  void
  print_strength_matrix(std::ostream &os = std::cout) const {
    os << "Strength matrix is: " << std::endl;
      for (unsigned int i = 0; i < strength_matrix.rows(); i++) {
          for (unsigned int j = 0; j < strength_matrix.cols(); j++) {
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
  print_C_i_set(const unsigned int i, std::ostream &os = std::cout) const {
    os << "Interpolating C-points of " << i << std::endl;
      for (auto c : map_point_to_interpolating_points.at(i)) {
        os << c << " - ";
      }
    os << std::endl;
  }

  void
  print_all_interpolating_points() const {
      for (auto i : F_points) {
        print_C_i_set(i);
      }
  }

  void
  press_to_continue(std::ostream &os = std::cout, std::istream &is = std::cin) const {
    os << "Press enter to continue: ";
    is.ignore();
  }

  //---------------------------------------------------------------------------------
  // PUBLIC MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
public:
  // Apply multigrid starting from physical domain
  AlgebraicMultigrid(const SparseMatrix &A, const std::vector<double> &b) : A(A), b(b) {
    std::cout << "Setting up AMG:" << std::endl;
    std::cout << "Matrix size: " << A.rows() << "x" << A.cols() << std::endl;

    construct_strength_sets();

      if (DEBUG) {
        print_strength_sets();
    }

    C_F_splitting();

      if (DEBUG) {
        print_C_points();
        print_F_points();
    }

    std::cout << "SETUP DONE!" << std::endl;
  }

  void
  solve(std::vector<double> &x) {
    std::cout << "Solving AMG:" << std::endl;
    unsigned int max_it   = 1000;
    int          tot_iter = 0;
    double       tol      = 1.e-6;
    double       normb;

    Jacobi jac(A, b, tol, nu);
    normb = Jacobi::norm(b);
      if (normb == 0.0) {
        normb = 1.0;
    }

    int flag = jac.solve(x);
    std::cout << "\tJAC1: flag " << flag << std::endl;

    tot_iter += jac.get_iter();
    std::cout << "\tJAC1: n_iter " << jac.get_iter() << std::endl;
    std::cout << "\tJAC1: tol_achieved " << jac.get_tol_achieved() << std::endl;

    // r_h = b - A*x_k - compute residual at step k on fine grid
    std::pair<std::vector<double>, bool> Ax_pair = A.mul(x);
    if (!Ax_pair.second)
      return;
    std::vector<double> r_h = jac.subvec(b, Ax_pair.first);

    auto r_2h = coarsen_residual(r_h);

    SparseMatrix A_2h = coarsen_matrix(A);

    // AlgebraicMultigrid  amg(D, A_2h, r_2h);
    std::vector<double> e_2h(r_2h.size(), 1.0);

      if (A.rows() > 3 && A.cols() > 3) {
        AlgebraicMultigrid inner_AMG(A_2h, r_2h);

        inner_AMG.solve(e_2h);
      } else {
        Jacobi coarse_Jac1(A_2h, r_2h, tol, max_it);

        flag = coarse_Jac1.solve(e_2h);
        std::cout << "\tCOARSEJAC1: flag " << flag << std::endl;

        tot_iter += coarse_Jac1.get_iter();
        std::cout << "\tCOARSEJAC1: n_iter " << coarse_Jac1.get_iter() << std::endl;
        std::cout << "\tCOARSEJAC1: tol_achieved " << coarse_Jac1.get_tol_achieved()
                  << std::endl;
        std::cout << "\ttot amg iter: " << tot_iter << std::endl;
      }

    auto e_h = interpolate_error(e_2h);

    jac.addvec_inplace(x, e_h);

    flag = jac.solve(x);
    std::cout << "\tJAC2: flag " << flag << std::endl;

    tot_iter += jac.get_iter();
    std::cout << "\tJAC2: n_iter " << jac.get_iter() << std::endl;
    std::cout << "\tJAC2: tol_achieved " << jac.get_tol_achieved() << std::endl;
    std::cout << "\ttot amg iter: " << tot_iter << std::endl;

    std::cout << "SOLVING: DONE!" << std::endl;
  }
};

#endif