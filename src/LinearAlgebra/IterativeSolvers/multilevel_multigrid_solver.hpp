#ifndef _MULTILEVEL_MULTIGRID_SOLVER_H__
#define _MULTILEVEL_MULTIGRID_SOLVER_H__

#include "multigrid_solver.hpp"

class MultilevelMultigridSolver : public MultigridSolver {
  //---------------------------------------------------------------------------------
  // PUBLIC MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
public:
  MultilevelMultigridSolver(const SparseMatrix   &A,
                            const Vector<double> &b,
                            const unsigned int   &pre_nu,
                            const unsigned int   &post_nu,
                            const double         &tol,
                            const unsigned int   &max_iter,
                            const unsigned int   &n_levels) :
    MultigridSolver(A, b, pre_nu, post_nu, tol, max_iter),
    n_levels(n_levels), system_matrices(n_levels), system_rhss(n_levels),
    restrictors(n_levels), interpolators(n_levels) {
    system_matrices[n_levels - 1] = std::make_unique<SparseMatrix>(A);
    system_rhss[n_levels - 1]     = std::make_unique<Vector<double>>(b);
  }

  virtual void
  build_all_restrictors() {
      for (int i = n_levels - 1; i > 0; i--) {
        build_restrictor(i);
      }
  }

  virtual void
  build_all_interpolators() {
      for (int i = n_levels - 1; i > 0; i--) {
        build_interpolator(i);
      }
  }

  virtual void
  coarsen_all_matrices() {
      for (int i = n_levels - 1; i > 0; i--) {
        build_restrictor(i);
        build_interpolator(i);

        coarsen_matrix(i);

        system_rhss[i - 1] = std::make_unique<Vector<double>>(restrictors[i]->rows());
      }
  }

  // multilevel multigrid solver
  virtual int
  solve(Vector<double> &x, const int &n_levels) = 0;
  //---------------------------------------------------------------------------------
  // PROTECTED MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
protected:
  unsigned int n_levels;

  Vector<std::unique_ptr<SparseMatrix>>   system_matrices;
  Vector<std::unique_ptr<Vector<double>>> system_rhss;
  Vector<std::unique_ptr<SparseMatrix>>   restrictors;
  Vector<std::unique_ptr<SparseMatrix>>   interpolators;

  virtual void
  build_restrictor(const unsigned int &lvl) = 0;

  virtual void
  build_interpolator(const unsigned int &lvl) = 0;

  virtual void
  coarsen_matrix(const unsigned int &lvl) = 0;
};

#endif