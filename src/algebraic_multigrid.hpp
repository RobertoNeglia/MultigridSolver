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
  const double pre_nu;
  const double post_nu;

  const double       tol;
  const unsigned int max_iter;

  double       tol_achieved;
  unsigned int n_iter;
  unsigned int jac_tot_iter;

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

  std::set<unsigned int> C_points;
  std::set<unsigned int> F_points;

  // Associates to each point i the c-points that are used to interpolate i
  std::map<unsigned int, std::set<unsigned int>> map_point_to_interpolating_points;

  template <typename T>
  unsigned int
  get_pos_in_set(std::set<T> set, T val) {
    unsigned int i = 0;
      for (auto it = set.begin(); it != set.end(); it++) {
        if (*it == val)
          return i;
        i++;
      }
    return set.size();
  }

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
    std::cout << "STRENGTH SETS: CONSTRUCTION..." << std::endl;
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

    max_lambda = 0;
      for (auto [i, s_i] : map_point_to_strength_transpose_set) {
        map_point_to_number_of_influences[i] = s_i.size();
        if (s_i.size() > max_lambda)
          max_lambda = s_i.size();
      }

    std::cout << "STRENGTH SETS: DONE!" << std::endl;
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
    std::cout << "CF SPLITTING: START..." << std::endl;
    const unsigned int n_elements = A.rows();
    // std::cout << "n_elements: " << n_elements << std::endl;
    // std::cout << "max_lambda: " << max_lambda << std::endl;
    // press_to_continue();
      while (C_points.size() + F_points.size() < n_elements) {
          for (auto [i, lambda_i] : map_point_to_number_of_influences) {
              if (lambda_i >= max_lambda && not_C_point(i) && not_F_point(i)) {
                // C point found
                // std::cout << "C point found! " << std::endl;
                // press_to_continue();
                make_C_point(i);
            }
          }
        if (max_lambda > 0)
          max_lambda--;
      }
    std::cout << "CF SPLITTING: DONE!" << std::endl;
  }

  std::vector<double>
  coarsen_residual(const std::vector<double> &r) {
    std::vector<double> rc(C_points.size());

    auto c_it = C_points.begin();

      for (unsigned int i = 0; i < rc.size(); i++) {
        rc[i] = r[*c_it];
        c_it++;
      }

    // double val;
    // for (unsigned int i = 0; i < rc.size(); i++) {
    //     if (i == 0) {
    //       val = (r[*c_it] + r[*c_it + 1]) / 2;
    //     } else if (i == rc.size() - 1) {
    //       val = (r[*c_it - 1] + r[*c_it]) / 2;
    //     } else {
    //       val = (r[*c_it - 1] + r[*c_it] + r[*c_it + 1]) / 3;
    //     }
    //   rc[i] = val;
    //   c_it++;
    // }

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

    auto j = e.begin();
      for (unsigned int i : C_points) {
        ef[i] = *j;
      }

    double       val = 0.0;
    unsigned int cnt = 0;
      for (unsigned int i : F_points) {
          for (unsigned int j : map_point_to_strength_set[i]) {
              if (C_points.contains(j)) {
                cnt++;
                val += e[get_pos_in_set(C_points, j)];
            }
          }
        val /= cnt;
        ef[i] = val;
      }

    // if (*(C_points.begin()) == 0) {
    //   std::cout << "PIPPOOOOOO" << std::endl;
    //     for (unsigned int i = 0; i < e.size() - 1; i++) {
    //       ef[2 * i]     = e[i];
    //       ef[2 * i + 1] = (e[i] + e[i + 1]) / 2;
    //     }
    //   ef[ef.size() - 1] = e[e.size() - 1];
    // } else {
    //   ef[0] = e[0];
    //   ef[1] = e[0];
    //     for (unsigned int i = 1; i < e.size(); i++) {
    //       ef[2 * i] = (e[i - 1] + e[i]) / 2;
    //       if (2 * i + 1 < ef.size() - 1)
    //         ef[2 * i + 1] = e[i];
    //     }
    // }

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

  void
  print_vector(std::vector<double> v) const {
      for (auto i : v) {
        std::cout << i << " - ";
      }
    std::cout << std::endl;
  }

  //---------------------------------------------------------------------------------
  // PUBLIC MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
public:
  // Apply multigrid starting from physical domain
  AlgebraicMultigrid(const SparseMatrix        &A,
                     const std::vector<double> &b,
                     const unsigned int         pre_nu,
                     const unsigned int         post_nu,
                     const double               tol,
                     const unsigned int         max_iter) :
    pre_nu(pre_nu),
    post_nu(post_nu), tol(tol), max_iter(max_iter), A(A), b(b) {
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

  int
  solve(std::vector<double> &x /*, int n_levels , std::vector<double> exac_sol*/) {
    if (DEBUG)
      std::cout << "Solving AMG:" << std::endl;
    int          tot_iter            = 0;
    double       inner_solver_tol    = 9 * 1.e-1;
    unsigned int inner_solver_max_it = 1000;
    double       normb;
    double       resid;

    normb = Jacobi::norm(b);
      if (normb == 0.0) {
        normb = 1.0;
    }

      for (unsigned int i = 0; i < max_iter; i++) {
        // pre-smoothing
        Jacobi pre_smoother(A, b, tol, pre_nu);
        int    flag = pre_smoother.solve(x);
        tot_iter += pre_smoother.get_iter();
          if (DEBUG) {
            std::cout << "\t PRE_SMOOTHER: flag " << flag << std::endl;
            std::cout << "\t PRE_SMOOTHER: n_iter " << pre_smoother.get_iter() << std::endl;
            std::cout << "\t PRE_SMOOTHER: tol_achieved " << pre_smoother.get_tol_achieved()
                      << std::endl;
        }

        // r_h = b - A*x_k - compute residual at step k on fine grid
        std::pair<std::vector<double>, bool> Ax_pair = A.mul(x);
        if (!Ax_pair.second)
          return -1;

        std::vector<double> r_h = Jacobi::subvec(b, Ax_pair.first);
        // print_vector(r_h);
        std::vector<double> r_2h = coarsen_residual(r_h);
        // print_vector(r_2h);

        // press_to_continue();

        SparseMatrix        A_2h = coarsen_matrix(A);
        std::vector<double> e_2h(r_2h.size(), 1.0);

        Jacobi coarse_smoother(A_2h, r_2h, inner_solver_tol, inner_solver_max_it);
        flag = coarse_smoother.solve(e_2h);
        tot_iter += coarse_smoother.get_iter();
          if (DEBUG) {
            std::cout << "\t COARSE_SOLVER: flag " << flag << std::endl;
            std::cout << "\t COARSE_SOLVER: n_iter " << coarse_smoother.get_iter()
                      << std::endl;
            std::cout << "\t COARSE_SOLVER: tol_achieved "
                      << coarse_smoother.get_tol_achieved() << std::endl;
        }

        // // multilevel
        //   if (n_levels > 0) {
        //     AlgebraicMultigrid inner_AMG(A_2h, r_2h, pre_nu, post_nu, tol, max_iter);
        //     inner_AMG.solve(e_2h, n_levels - 1);
        //   } else {
        //     Jacobi coarse_smoother(A_2h, r_2h, inner_solver_tol, inner_solver_max_it);
        //     flag = coarse_smoother.solve(e_2h);
        //     tot_iter += coarse_smoother.get_iter();
        //       if (DEBUG) {
        //         std::cout << "\t COARSE_SOLVER: flag " << flag << std::endl;
        //         std::cout << "\t COARSE_SOLVER: n_iter " << coarse_smoother.get_iter()
        //                   << std::endl;
        //         std::cout << "\t COARSE_SOLVER: tol_achieved "
        //                   << coarse_smoother.get_tol_achieved() << std::endl;
        //     }
        //   }

        std::vector<double> e_h = interpolate_error(e_2h);
        // print_vector(e_h);
        // press_to_continue();

        Jacobi::addvec_inplace(x, e_h);

        // std::vector<double> e = Jacobi::subvec(exac_sol, x);
        // print_vector(e);

        // post-smoothing
        Jacobi post_smoother(A, b, tol, post_nu);
        flag = post_smoother.solve(x);
        tot_iter += post_smoother.get_iter();
          if (DEBUG) {
            std::cout << "\t POST_SMOOTHER: flag " << flag << std::endl;
            std::cout << "\t POST_SMOOTHER: n_iter " << post_smoother.get_iter() << std::endl;
            std::cout << "\t POST_SMOOTHER: tol_achieved " << post_smoother.get_tol_achieved()
                      << std::endl;
        }

        std::pair<std::vector<double>, bool> Ax_k_pair = A.mul(x);
        if (!Ax_k_pair.second)
          return -1;
        std::vector<double> r_k = Jacobi::subvec(b, Ax_k_pair.first);

          if ((resid = Jacobi::norm(r_k) / normb) <= tol) {
            tol_achieved = resid;
            n_iter       = i;
            jac_tot_iter = tot_iter;
            std::cout << "TOLERANCE ACHIEVED" << std::endl;
            std::cout << "SOLVING: DONE!" << std::endl;
            return 1;
        }
      }

    tol_achieved = resid;
    n_iter       = max_iter;
    jac_tot_iter = tot_iter;
    std::cout << "ITERATION LIMIT REACHED" << std::endl;
    std::cout << "SOLVING: DONE!" << std::endl;
    return 2;
  }

  unsigned int
  get_iter() const {
    return n_iter;
  }

  double
  get_tol_achieved() const {
    return tol_achieved;
  }

  unsigned int
  get_jac_tot_iter() const {
    return jac_tot_iter;
  }
};

#endif