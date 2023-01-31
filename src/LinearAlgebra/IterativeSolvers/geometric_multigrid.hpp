#ifndef _GEOMETRIC_MULTIGRID_H__
#define _GEOMETRIC_MULTIGRID_H__

#include "../sparse_matrix.hpp"
#include "jacobi_r.hpp"
#include "multigrid_solver.hpp"

#define DEBUG_SOLVING 0

class GeometricMultigrid : public MultigridSolver {
  //---------------------------------------------------------------------------------
  // PUBLIC MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
public:
  using MultigridSolver::MultigridSolver;

  void
  setup() override {
    std::cout << "==========================================================" << std::endl;
    std::cout << "Setting up GMG: " << std::endl;
    std::cout << "Matrix size: " << A.rows() << "x" << A.cols() << std::endl;

    build_restrictor();
    build_interpolator();

    coarsen_matrix(A);

    /* DEBUG
      if (DEBUG_SOLVING) {
        A_2h.print_matrix();
        press_to_continue();
    }

    if (DEBUG_INTERPOLATOR) {
            std::cout << "Restrictor size: " << restrictor.rows() << "x" << restrictor.cols()
                      << std::endl;
            // restrictor.print_matrix();
            std::cout << "Interpolator size: " << interpolator.rows() << "x" <<
       interpolator.cols()
                      << std::endl;
            // interpolator.print_matrix();
    }
    */

    std::cout << "Setup DONE!" << std::endl << std::endl;
  }

  int
  solve(Vector<double> &x) override {
    std::cout << "==========================================================" << std::endl;
    std::cout << "Solving GMG:" << std::endl;

    unsigned int tot_iter               = 0;
    double       coarse_smoother_tol    = 1.e-6;
    unsigned int coarse_smoother_max_it = 100;
    double       normb;
    double       resid;

    Vector<double> b_k(x.size());
    Vector<double> r_h(x.size());
    Vector<double> r_2h(restrictor.rows());
    Vector<double> e_2h(restrictor.rows());
    Vector<double> e_h(x.size());

    normb = norm(b);
      if (normb == 0.0) {
        normb = 1.0;
    }

    pre_smoother = std::make_unique<Jacobi>(A, b, tol, pre_nu);
    coarse_smoother =
      std::make_unique<Jacobi>(A_2h, r_2h, coarse_smoother_tol, coarse_smoother_max_it);
    post_smoother = std::make_unique<Jacobi>(A, b, tol, post_nu);

      for (unsigned int i = 0; i < max_iter; i++) {
        // pre-smoothing
        int flag = pre_smoother->solve(x);
        tot_iter += pre_smoother->get_iter();

        /* DEBUG
          if (DEBUG_SOLVING) {
            std::cout << "  PRE_SMOOTHER: flag " << flag << std::endl;
            std::cout << "  PRE_SMOOTHER: n_iter " << pre_smoother->get_iter() << std::endl;
            std::cout << "  PRE_SMOOTHER: tol_achieved " << pre_smoother->get_tol_achieved()
                      << std::endl;
        }
        */

        // compute residual at step k on fine grid
        A.mul(b_k, x);       // A * x_k
        subvec(r_h, b, b_k); // b - A * x_k

        // bring residual to coraser grid
        restrictor.mul(r_2h, r_h);

        // coarse grid solver inital guess
        fill(e_2h, 0.0);

        // solve the error equation on the coarse grid
        coarse_smoother->change_rhs(r_2h);
        flag = coarse_smoother->solve(e_2h);
        tot_iter += coarse_smoother->get_iter();

        /* DEBUG
          if (DEBUG_SOLVING) {
            std::cout << "  COARSE_SOLVER: flag " << flag << std::endl;
            std::cout << "  COARSE_SOLVER: n_iter " << coarse_smoother->get_iter()
                      << std::endl;
            std::cout << "  COARSE_SOLVER: tol_achieved "
                      << coarse_smoother->get_tol_achieved() << std::endl;
        }
        */

        // bring error back to the fine grid
        interpolator.mul(e_h, e_2h);

        // update approximate solution with the error computed on the fine grid
        addvec_inplace(x, e_h);

        // post-smoothing
        flag = post_smoother->solve(x);
        tot_iter += post_smoother->get_iter();

        /* DEBUG
          if (DEBUG_SOLVING) {
            std::cout << "  POST_SMOOTHER: flag " << flag << std::endl;
            std::cout << "  POST_SMOOTHER: n_iter " << post_smoother->get_iter() <<
            std::endl; std::cout << "  POST_SMOOTHER: tol_achieved " <<
            post_smoother->get_tol_achieved()
                      << std::endl;
        }
        */

        subvec(r_h, b, A.mul(x));

          if ((resid = norm(r_h) / normb) <= tol) {
            tol_achieved      = resid;
            n_iter            = i;
            smoother_tot_iter = tot_iter;
            pre_smoother.reset();
            coarse_smoother.reset();
            post_smoother.reset();
            std::cout << "  TOLERANCE ACHIEVED" << std::endl;
            std::cout << "SOLVING: DONE!" << std::endl << std::endl;
            return 1;
        }
      }
    tol_achieved      = resid;
    n_iter            = max_iter;
    smoother_tot_iter = tot_iter;
    pre_smoother.reset();
    coarse_smoother.reset();
    post_smoother.reset();
    std::cout << "  ITERATION LIMIT REACHED" << std::endl;
    std::cout << "SOLVING: DONE!" << std::endl << std::endl;
    return 2;
  }

  //---------------------------------------------------------------------------------
  // PRIVATE MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
private:
  void
  build_restrictor() override {
    std::cout << "  Building the restrictor..." << std::endl;
    unsigned int N = std::sqrt(A.rows());
    unsigned int M = std::floor((N - 1) / 2);
    unsigned int m = M * M;
    unsigned int n = A.cols();

    restrictor.initialize(m, n);
      for (unsigned int il = 0; il < M; il++) {
        unsigned int i = 1 + il * 2;
          for (unsigned int jl = 0; jl < M; jl++) {
            unsigned int j   = 1 + jl * 2;
            unsigned int ijl = il + jl * M;

            restrictor.insert_coeff(1.0 / 4.0, ijl, i + j * N);
            restrictor.insert_coeff(1.0 / 8.0, ijl, (i - 1) + j * N);
            restrictor.insert_coeff(1.0 / 8.0, ijl, (i + 1) + j * N);
            restrictor.insert_coeff(1.0 / 8.0, ijl, i + (j - 1) * N);
            restrictor.insert_coeff(1.0 / 8.0, ijl, i + (j + 1) * N);
            restrictor.insert_coeff(1.0 / 16.0, ijl, (i - 1) + (j - 1) * N);
            restrictor.insert_coeff(1.0 / 16.0, ijl, (i - 1) + (j + 1) * N);
            restrictor.insert_coeff(1.0 / 16.0, ijl, (i + 1) + (j - 1) * N);
            restrictor.insert_coeff(1.0 / 16.0, ijl, (i + 1) + (j + 1) * N);
          }
      }

    std::cout << "  DONE!" << std::endl;
  }

  void
  build_interpolator() override {
    std::cout << "  Building the interpolator..." << std::endl;

    interpolator = restrictor.transpose();

    std::cout << "  DONE!" << std::endl;
  }
};

#endif