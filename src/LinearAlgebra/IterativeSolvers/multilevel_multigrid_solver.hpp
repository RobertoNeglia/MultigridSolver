#ifndef _MULTILEVEL_MULTIGRID_SOLVER_H__
#define _MULTILEVEL_MULTIGRID_SOLVER_H__

#include "multigrid_solver.hpp"

class MultilevelMultigridSolver : public MultigridSolver {
  //---------------------------------------------------------------------------------
  // PUBLIC MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
public:
  MultilevelMultigridSolver(SparseMatrix      &A,
                            Vector<double>    &b,
                            const unsigned int pre_nu,
                            const unsigned int post_nu,
                            const double       tol,
                            const unsigned int max_iter,
                            const unsigned int n_levels) :
    MultigridSolver(A, b, pre_nu, post_nu, tol, max_iter),
    n_levels(n_levels) {
    system_matrices.reserve(n_levels);
    system_rhss.reserve(n_levels);
    restrictors.reserve(n_levels);
    interpolators.reserve(n_levels);
    system_matrices[n_levels - 1] = &A;
    system_rhss[n_levels - 1]     = &b;
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

        system_rhss[i - 1] = new Vector<double>(restrictors[i]->rows());
      }
  }

  virtual void
  free_space() {
      for (auto i : system_matrices) {
        free(i);
      }

      for (auto i : system_rhss) {
        free(i);
      }

      for (auto i : restrictors) {
        free(i);
      }

      for (auto i : interpolators) {
        free(i);
      }
  }

  // multilevel multigrid solver
  virtual int
  solve(Vector<double> &x, const int n_levels) = 0;
  //---------------------------------------------------------------------------------
  // PROTECTED MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
protected:
  unsigned int n_levels;

  Vector<SparseMatrix *>   system_matrices;
  Vector<Vector<double> *> system_rhss;
  Vector<SparseMatrix *>   restrictors;
  Vector<SparseMatrix *>   interpolators;

  virtual void
  build_restrictor(const unsigned int lvl) = 0;

  virtual void
  build_interpolator(const unsigned int lvl) = 0;

  virtual void
  coarsen_matrix(const unsigned int lvl) = 0;
};

#endif