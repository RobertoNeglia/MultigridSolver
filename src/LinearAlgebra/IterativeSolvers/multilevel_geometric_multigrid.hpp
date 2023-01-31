#ifndef _MULTILEVEL_GEOMETRIC_MULTIGRID_H__
#define _MULTILEVEL_GEOMETRIC_MULTIGRID_H__

#include "jacobi_r.hpp"
#include "multilevel_multigrid_solver.hpp"

#define DEBUG_SOLVING 0

class MultilevelGeometricMultigrid : public MultilevelMultigridSolver {
  //---------------------------------------------------------------------------------
  // PUBLIC MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
public:
  using MultilevelMultigridSolver::MultilevelMultigridSolver;

  void
  setup() override {
    std::cout << "==========================================================" << std::endl;
    std::cout << "Setting up GMG: " << std::endl;
    std::cout << "Matrix size: " << A.rows() << "x" << A.cols() << std::endl;

    coarsen_all_matrices();

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

  virtual int
  solve(Vector<double> &x) override {
    UNUSED(x);
    return 0;
  }

  virtual int
  solve(Vector<double> &x, const int curr_lvl) override {
    const int lvl = curr_lvl - 1;

    std::cout << "==========================================================" << std::endl;
    std::cout << std::string((n_levels - lvl) * 2, ' ') << "Solving GMG:" << std::endl;
    std::cout << std::string((n_levels - lvl) * 2, ' ') << "LEVEL : " << (curr_lvl)
              << std::endl;

    unsigned int tot_iter               = 0;
    double       coarse_smoother_tol    = 1.e-6;
    unsigned int coarse_smoother_max_it = 50;
    double       normb;
    double       resid;

    Vector<double> b_k(x.size());
    Vector<double> r_h(x.size());
    Vector<double> e_2h(restrictors[lvl]->rows());
    Vector<double> e_h(x.size());

    normb = norm(b);
      if (normb == 0.0) {
        normb = 1.0;
    }

      for (unsigned int i = 0; i < max_iter; i++) {
        pre_smoother = std::make_unique<Jacobi>(*system_matrices[lvl], // A^i
                                                *system_rhss[lvl],     // b^i
                                                tol,
                                                pre_nu);
        // pre-smoothing
        int flag = pre_smoother->solve(x);
        tot_iter += pre_smoother->get_iter();

        // compute residual at step k on fine grid
        system_matrices[lvl]->mul(b_k, x);   // A^i * x^i_k
        subvec(r_h, *system_rhss[lvl], b_k); // b^i - A^i * x^i_k = r^i

        // bring residual to coraser grid
        restrictors[lvl]              // I^i-1_i - restrictor from level i to level i-1
          ->mul(system_rhss[lvl - 1], // r^(i-1)
                r_h);                 // r^i

        // coarse grid solver inital guess
        fill(e_2h, 0.0);

          if (lvl >= 2 && system_matrices[lvl]->rows() >= 9) {
            std::cout << std::string((n_levels - lvl + 1) * 2, ' ')
                      << "CALLING ANOTHER LEVEL OF MULTIGRID..." << std::endl;
            solve(e_2h, curr_lvl - 1); // MG(e^i-1, b^i-1)
            std::cout << std::string((n_levels - lvl + 1) * 2, ' ') << "BACK TO GRID LEVEL "
                      << curr_lvl << std::endl;
          } else {
            // solve the error equation on the coarser grid
            if (!i)
              std::cout << std::string((n_levels - lvl + 1) * 2, ' ')
                        << "SOLVING THE ERROR EQUATION ON THE COARSEST GRID..." << std::endl;
            coarse_smoother = std::make_unique<Jacobi>(*system_matrices[lvl - 1],
                                                       *system_rhss[lvl - 1],
                                                       coarse_smoother_tol,
                                                       coarse_smoother_max_it);
            flag            = coarse_smoother->solve(e_2h);
            tot_iter += coarse_smoother->get_iter();
          }

        // bring error back to the fine grid
        interpolators[lvl] // I^i_i-1 - interpolator from level i-1 to level i
          ->mul(e_h,       // e^i
                e_2h);     // e^i-1

        // update approximate solution with the error computed on the fine grid
        addvec_inplace(x, e_h);

        post_smoother = std::make_unique<Jacobi>(*system_matrices[lvl], // A^i
                                                 *system_rhss[lvl],     // b^i
                                                 tol,
                                                 post_nu);

        // post-smoothing
        flag = post_smoother->solve(x);
        tot_iter += post_smoother->get_iter();

        subvec(r_h, *system_rhss[lvl], system_matrices[lvl]->mul(x));

          if ((resid = norm(r_h) / normb) <= tol) {
            tol_achieved      = resid;
            n_iter            = i;
            smoother_tot_iter = tot_iter;
            pre_smoother.reset();
            coarse_smoother.reset();
            post_smoother.reset();
            std::cout << std::string((n_levels - lvl + 1) * 2, ' ') << "TOLERANCE ACHIEVED"
                      << std::endl;
            std::cout << std::string((n_levels - lvl) * 2, ' ') << "SOLVING: DONE!"
                      << std::endl
                      << std::endl;
            return 1;
        }
      }
    tol_achieved      = resid;
    n_iter            = max_iter;
    smoother_tot_iter = tot_iter;
    pre_smoother.reset();
    coarse_smoother.reset();
    post_smoother.reset();

    std::cout << std::string((n_levels - lvl + 1) * 2, ' ') << "ITERATION LIMIT REACHED"
              << std::endl;
    std::cout << std::string((n_levels - lvl) * 2, ' ') << "SOLVING: DONE!" << std::endl
              << std::endl;
    return 2;
  }

  //---------------------------------------------------------------------------------
  // PRIVATE MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
private:
  void
  build_restrictor(const unsigned int lvl) override {
    std::cout << "  Building the restrictor of level " << lvl << "..." << std::endl;

    unsigned int N = std::sqrt(system_matrices[lvl]->rows());
    unsigned int M = std::floor((N - 1) / 2);
    unsigned int m = M * M;
    unsigned int n = system_matrices[lvl]->cols();

    SparseMatrix *restrictor = new SparseMatrix;

    restrictor->initialize(m, n);
      for (unsigned int il = 0; il < M; il++) {
        unsigned int i = 1 + il * 2;
          for (unsigned int jl = 0; jl < M; jl++) {
            unsigned int j   = 1 + jl * 2;
            unsigned int ijl = il + jl * M;

            restrictor->insert_coeff(1.0 / 4.0, ijl, i + j * N);
            restrictor->insert_coeff(1.0 / 8.0, ijl, (i - 1) + j * N);
            restrictor->insert_coeff(1.0 / 8.0, ijl, (i + 1) + j * N);
            restrictor->insert_coeff(1.0 / 8.0, ijl, i + (j - 1) * N);
            restrictor->insert_coeff(1.0 / 8.0, ijl, i + (j + 1) * N);
            restrictor->insert_coeff(1.0 / 16.0, ijl, (i - 1) + (j - 1) * N);
            restrictor->insert_coeff(1.0 / 16.0, ijl, (i - 1) + (j + 1) * N);
            restrictor->insert_coeff(1.0 / 16.0, ijl, (i + 1) + (j - 1) * N);
            restrictor->insert_coeff(1.0 / 16.0, ijl, (i + 1) + (j + 1) * N);
          }
      }

    restrictors[lvl] = restrictor;

    std::cout << "  DONE!" << std::endl;
  }

  void
  build_restrictor() override {}

  void
  build_interpolator() override {}

  void
  build_interpolator(const unsigned int lvl) override {
    std::cout << "  Building the interpolator of level " << lvl << "..." << std::endl;

    interpolators[lvl] = &(restrictors[lvl]->transpose());

    std::cout << "  DONE!" << std::endl;
  }

  void
  coarsen_matrix(const unsigned int lvl) {
    std::cout << "  Coarsening the system matrix to level " << (lvl - 1) << "..." << std::endl;

    system_matrices[lvl - 1]          //
      = restrictors[lvl]              //
          ->mul(system_matrices[lvl]) //
          ->mul(interpolators[lvl]);

    std::cout << "  Coarsening done" << std::endl
              << "  Coarsen matrix size: " << system_matrices[lvl - 1]->rows() << "x"
              << system_matrices[lvl - 1]->cols() << std::endl;
  }
};

#endif