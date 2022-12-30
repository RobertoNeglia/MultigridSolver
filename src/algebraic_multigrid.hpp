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
  // Alpha parameter that determines strong connection
  const double alpha = 1.0 / 4.0;

  // number of iterations before & after multigrid correction
  const double nu = 1000;

  // Domain of the problem
  Domain2D D;

  // Matrix A of the linear system
  SparseMatrix A;

  // System rhs
  std::vector<double> b;

  // Strength set - to each point i is associated the set of points j that influence i
  //                  SET OF POINTS j THAT INFLUENCE i
  std::map<int, std::set<unsigned int>> map_point_to_strength_set;

  // Strength transpose set - to each point i is associated the set of points j that are
  // influenced by i
  //                  SET OF POINTS j THAT ARE INFLUENCED BY i
  std::map<int, std::set<unsigned int>> map_point_to_strength_transpose_set;

  // Associates to each point i the number of points that are influenced by i
  std::map<unsigned int, unsigned int> map_point_to_number_of_influences;

  int max_lambda;

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
      }
  }

  void
  construct_strength_matrix() {
    unsigned int max = 0;
      for (auto [i, lambda_i] : map_point_to_number_of_influences) {
        std::pair<unsigned int, unsigned int> coordinates = D.onedim2twodim(i);
        strength_matrix.insert_coeff(lambda_i, coordinates.first, coordinates.second);
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
        int &lambda = strength_matrix.coeff_ref(i, j);
        // reinterpret_cast<unsigned int &>(strength_matrix.coeffRef(i, j));
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
  }

  void
  make_F_point(const unsigned int i, const unsigned int j, const unsigned int c_point) {
      if (inside_boundary(i, j) && not_C_point(i, j)) {
        const unsigned int f_point = D.twodim2onedim(i, j);
        map_point_to_interpolating_points[f_point].insert(c_point);
          if (not_F_point(i, j)) {
            F_points.insert(f_point);
            update_neighbours_strength(i, j);
        }
    }
  }

  void
  neighbours_to_F_points(const unsigned int i, const unsigned int j) {
    const unsigned int c_point = D.twodim2onedim(i, j);
    make_F_point(i - 1, j, c_point);
    make_F_point(i, j - 1, c_point);
    make_F_point(i + 1, j, c_point);
    make_F_point(i, j + 1, c_point);
  }

  void
  make_C_point(const unsigned int i, const unsigned int j) {
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
  C_F_splitting() {
    const unsigned int n_elements = D.rows() * D.cols();

      while (C_points.size() + F_points.size() < n_elements) {
          for (unsigned int i = 0; i < strength_matrix.rows(); i++) {
              for (unsigned int j = 0; j < strength_matrix.cols(); j++) {
                int lambda = strength_matrix.coeff(i, j);
                  if (lambda >= max_lambda && not_C_point(i, j) && not_F_point(i, j)) {
                    make_C_point(i, j);
                }
              }
          }
        max_lambda--;
      }
  }

  void
  coarsen_domain() {
      for (auto i : C_points) {
        std::pair<unsigned int, unsigned int> ij = D.onedim2twodim(i);
        coarsened_domain.insert_coeff(i, ij.first, ij.second);
      }

    coarsened_domain.print();
  }

  std::vector<double>
  coarsen_residual(const std::vector<double> &r) {
    std::vector<double> rc(C_points.size());
      for (unsigned int i = 0; i < rc.size(); i++) {
        unsigned int j = 2 * i;
        if (j < r.size() - 1)
          rc[i] = 0.5 * (r[j] + r[j + 1]);
        else if (j == r.size() - 1)
          rc[i] = r[j];
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
            double       new_coeff;
            // Boundary management
              if (ii < A.rows() - 1 && jj < A.cols() - 1) {
                new_coeff =
                  0.25 *
                  (A.coeff(2 * i, 2 * j).first + A.coeff(2 * i + 1, 2 * j).first +
                   A.coeff(2 * i, 2 * j + 1).first + A.coeff(2 * i + 1, 2 * j + 1).first);
              } else {
                  if (ii == A.rows() - 1) {
                    new_coeff =
                      0.5 * (A.coeff(2 * i, 2 * j).first + A.coeff(2 * i, 2 * j + 1).first);
                }
                  if (jj == A.cols() - 1) {
                    new_coeff =
                      0.5 * (A.coeff(2 * i, 2 * j).first + A.coeff(2 * i + 1, 2 * j).first);
                }
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
  AlgebraicMultigrid(const Domain2D &D, const SparseMatrix &A, const std::vector<double> &b) :
    D(D), A(A), b(b), strength_matrix(D.rows(), D.cols()) {
    std::cout << "Setting up AMG:" << std::endl;
    std::cout << "Matrix size: " << A.rows() << "x" << A.cols() << std::endl;
    construct_strength_sets();
    construct_strength_matrix();
      if (DEBUG) {
        print_strength_sets();
        print_strength_matrix();
    }

    C_F_splitting();
    // coarsen_domain();

      if (DEBUG) {
        print_C_points();
        print_F_points();
    }

    SparseMatrix Ac = coarsen_matrix(A);
    // Ac.print_matrix();
    std::cout << "DONE!" << std::endl;
  }

  AlgebraicMultigrid(const std::set<unsigned int> &c_points,
                     const SparseMatrix           &A,
                     const std::vector<double>    &b) :
    A(A),
    b(b) {
    std::cout << "Setting up AMG:" << std::endl;
    std::cout << "Matrix size: " << A.rows() << "x" << A.cols() << std::endl;
    construct_strength_sets();
    construct_strength_matrix();
      if (DEBUG) {
        print_strength_sets();
        print_strength_matrix();
    }

    C_F_splitting();

      if (DEBUG) {
        print_C_points();
        print_F_points();
    }

    SparseMatrix Ac = coarsen_matrix(A);
    // Ac.print_matrix();
    std::cout << "DONE!" << std::endl;
  }

  void
  solve(std::vector<double> &x) {
    std::cout << "Solving AMG:" << std::endl;
    unsigned int max_it = 2000;
    double       tol    = 1.e-6;
    Jacobi       jac(A, b, tol, nu);

    int flag = jac.solve(x);
    std::cout << "JAC1: flag " << flag << std::endl;

    int tot_iter = jac.get_iter();
    std::cout << "JAC1: n_iter " << jac.get_iter() << std::endl;
    std::cout << "JAC1: tol_achieved " << jac.get_tol_achieved() << std::endl;

    // r_h = b - A*x_k - compute residual at step k on fine grid
    std::pair<std::vector<double>, bool> Ax_pair = A.mul(x);
    if (!Ax_pair.second)
      return;
    std::vector<double> r_h = jac.subvec(b, Ax_pair.first);

    auto r_2h = coarsen_residual(r_h);

    SparseMatrix A_2h = coarsen_matrix(A);

    // AlgebraicMultigrid  amg(D, A_2h, r_2h);
    std::vector<double> e_2h(r_2h.size(), 1.0);

    Jacobi coarse_Jac1(A_2h, r_2h, tol, max_it);

    flag = coarse_Jac1.solve(e_2h);
    std::cout << "COARSEJAC1: flag " << flag << std::endl;

    tot_iter += coarse_Jac1.get_iter();
    std::cout << "COARSEJAC1: n_iter " << coarse_Jac1.get_iter() << std::endl;
    std::cout << "COARSEJAC1: tol_achieved " << coarse_Jac1.get_tol_achieved() << std::endl;
    std::cout << "tot amg iter: " << tot_iter << std::endl;

    auto e_h = interpolate_error(e_2h);

    jac.addvec_inplace(x, e_h);

    flag = jac.solve(x);
    std::cout << "JAC2: flag " << flag << std::endl;

    tot_iter += jac.get_iter();
    std::cout << "JAC2: n_iter " << jac.get_iter() << std::endl;
    std::cout << "JAC2: tol_achieved " << jac.get_tol_achieved() << std::endl;
    std::cout << "tot amg iter: " << tot_iter << std::endl;

    /*
        Ax_pair = A.mul(x);
        if (!Ax_pair.second)
          return;

        r_h  = jac.subvec(b, Ax_pair.first);
        r_2h = coarsen_residual(r_h);
        e_2h.resize(r_2h.size(), 1.0);

        Jacobi coarse_Jac2(A_2h, r_2h, tol, max_it);

        flag = coarse_Jac2.solve(e_2h);
        std::cout << "COARSEJAC2: flag " << flag << std::endl;

        tot_iter += coarse_Jac2.get_iter();
        std::cout << "COARSEJAC2: n_iter " << coarse_Jac2.get_iter() << std::endl;
        std::cout << "COARSEJAC2: tol_achieved " << coarse_Jac2.get_tol_achieved() <<
        std::endl; std::cout << "tot amg iter: " << tot_iter << std::endl;

        e_h = interpolate_error(e_2h);
        x   = jac.addvec_inplace(x, e_h);

        flag = jac.solve(x);
        std::cout << "JAC2: flag " << flag << std::endl;

        tot_iter += jac.get_iter();
        std::cout << "JAC2: n_iter " << jac.get_iter() << std::endl;
        std::cout << "JAC2: tol_achieved " << jac.get_tol_achieved() << std::endl;
        std::cout << "tot amg iter: " << tot_iter << std::endl;

        std::cout << "DONE!" << std::endl;
    */
  }
};

#endif