#ifndef __ALGEBRAIC_MULTIGRID_H_
#define __ALGEBRAIC_MULTIGRID_H_

#include <map>
#include <set>

#include "../sparse_matrix.hpp"
#include "jacobi_r.hpp"
#include "multigrid_solver.hpp"

#define DEBUG_CF_SPLITTING 0
#define DEBUG_SOLVING 0
#define DBL_MIN -100000000

class AlgebraicMultigrid : public MultigridSolver {
  //---------------------------------------------------------------------------------
  // PUBLIC MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
public:
  using MultigridSolver::MultigridSolver;

  void
  setup() override {
    std::cout << "Setting up AMG:" << std::endl;
    std::cout << "  Matrix size: " << A.rows() << "x" << A.cols() << std::endl;

    construct_strength_sets();

      if (DEBUG_CF_SPLITTING) {
        print_strength_sets();
    }

    C_F_splitting();

      if (DEBUG_CF_SPLITTING) {
        print_C_points();
        print_F_points();
    }

    build_interpolator();
    build_restrictor();

      if (DEBUG_INTERPOLATOR) {
        std::cout << "Restrictor size: " << restrictor.rows() << "x" << restrictor.cols()
                  << std::endl;
        // restrictor.print_matrix();
        std::cout << "Interpolator size: " << interpolator.rows() << "x" << interpolator.cols()
                  << std::endl;
        // interpolator.print_matrix();
        press_to_continue();
    }

    std::cout << "SETUP DONE!" << std::endl << std::endl;
  }

  int
  solve(std::vector<double> &x) override {
    std::cout << "Solving AMG:" << std::endl;
    int          tot_iter               = 0;
    double       coarse_smoother_tol    = 1.e-100;
    unsigned int coarse_smoother_max_it = 5;
    double       normb;
    double       resid;

    normb = norm(b);
      if (normb == 0.0) {
        normb = 1.0;
    }

    SparseMatrix A_2h = coarsen_matrix(A);
      if (DEBUG_SOLVING) {
        A_2h.print_matrix();
        press_to_continue();
    }

    pre_smoother  = std::make_unique<Jacobi>(A, b, tol, pre_nu);
    post_smoother = std::make_unique<Jacobi>(A, b, tol, post_nu);

      for (unsigned int i = 0; i < max_iter; i++) {
        // pre-smoothing
        int flag = pre_smoother->solve(x);
        tot_iter += pre_smoother->get_iter();
          if (DEBUG_SOLVING) {
            std::cout << "  PRE_SMOOTHER: flag " << flag << std::endl;
            std::cout << "  PRE_SMOOTHER: n_iter " << pre_smoother->get_iter() << std::endl;
            std::cout << "  PRE_SMOOTHER: tol_achieved " << pre_smoother->get_tol_achieved()
                      << std::endl;
        }

        // r_h = b - A*x_k - compute residual at step k on fine grid
        std::vector<double> r_h  = subvec(b, A.mul(x));
        std::vector<double> r_2h = restrictor.mul(r_h);

          if (DEBUG_SOLVING) {
            print_vector(r_h);
            print_vector(r_2h);
            press_to_continue();
        }

        std::vector<double> e_2h(r_2h.size(), 0.0);

        coarse_smoother =
          std::make_unique<Jacobi>(A_2h, r_2h, coarse_smoother_tol, coarse_smoother_max_it);
        flag = coarse_smoother->solve(e_2h);
        tot_iter += coarse_smoother->get_iter();
          if (DEBUG_SOLVING) {
            std::cout << "  COARSE_SOLVER: flag " << flag << std::endl;
            std::cout << "  COARSE_SOLVER: n_iter " << coarse_smoother->get_iter()
                      << std::endl;
            std::cout << "  COARSE_SOLVER: tol_achieved "
                      << coarse_smoother->get_tol_achieved() << std::endl;
        }

        std::vector<double> e_h = interpolator.mul(e_2h);

          if (DEBUG_SOLVING) {
            print_vector(e_2h);
            print_vector(e_h);
            press_to_continue();
        }

        addvec_inplace(x, e_h);

        // post-smoothing
        flag = post_smoother->solve(x);
        tot_iter += post_smoother->get_iter();
          if (DEBUG_SOLVING) {
            std::cout << "  POST_SMOOTHER: flag " << flag << std::endl;
            std::cout << "  POST_SMOOTHER: n_iter " << post_smoother->get_iter() << std::endl;
            std::cout << "  POST_SMOOTHER: tol_achieved " << post_smoother->get_tol_achieved()
                      << std::endl;
        }

        std::vector<double> r_k = subvec(b, A.mul(x));

          if ((resid = norm(r_k) / normb) <= tol) {
            tol_achieved      = resid;
            n_iter            = i;
            smoother_tot_iter = tot_iter;
            std::cout << "  TOLERANCE ACHIEVED" << std::endl;
            std::cout << "SOLVING: DONE!" << std::endl << std::endl;
            return 1;
        }
      }

    tol_achieved      = resid;
    n_iter            = max_iter;
    smoother_tot_iter = tot_iter;
    std::cout << "  ITERATION LIMIT REACHED" << std::endl;
    std::cout << "SOLVING: DONE!" << std::endl << std::endl;
    return 2;
  }

  virtual int
  solve(Vector<double> &x, unsigned int n_levels) override {
    UNUSED(x);
    UNUSED(n_levels);
    return 0;
  }

  //---------------------------------------------------------------------------------
  // PRIVATE MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
private:
  /**
   * Alpha parameter to determine strong connection
   */
  const double alpha = 1.0 / 4.0;

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
    std::cout << "  STRENGTH SETS: CONSTRUCTION..." << std::endl;
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

    std::cout << "  STRENGTH SETS: DONE!" << std::endl;
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
      if (DEBUG_CF_SPLITTING) {
        print_C_points();
        print_F_points();
    }
  }

  void
  C_F_splitting() {
    std::cout << "  CF SPLITTING: START..." << std::endl;
    const unsigned int n_elements = A.rows();
      while (C_points.size() + F_points.size() < n_elements) {
          for (auto [i, lambda_i] : map_point_to_number_of_influences) {
              if (lambda_i >= max_lambda && not_C_point(i) && not_F_point(i)) {
                make_C_point(i);
            }
          }
        if (max_lambda > 0)
          max_lambda--;
      }
    std::cout << "  CF SPLITTING: DONE!" << std::endl;
  }

  void
  build_interpolator() override {
    interpolator.initialize(C_points.size() + F_points.size(), C_points.size());

    unsigned int j = 0;
      for (auto i : C_points) {
        interpolator.insert_coeff(1.0, i, j);
        j++;
      }

      for (auto i : F_points) {
          for (unsigned int j : map_point_to_strength_set[i]) {
            interpolator.insert_coeff(1.0 / map_point_to_strength_set[i].size(),
                                      i,
                                      get_pos_in_set(C_points, j));
          }
      }
  }

  void
  build_restrictor() override {
    restrictor = interpolator.transpose();
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
  print_vector(std::vector<double> v) const {
      for (auto i : v) {
        std::cout << i << " - ";
      }
    std::cout << std::endl;
  }
};

#endif