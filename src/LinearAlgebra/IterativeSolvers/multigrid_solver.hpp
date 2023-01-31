#ifndef _MULTIGRID_SOLVER_H__
#define _MULTIGRID_SOLVER_H__

#include "iterative_solver.hpp"
#define DEBUG_INTERPOLATOR 0

#define UNUSED(expr) \
    do {             \
      (void)(expr);  \
  } while (0)

class MultigridSolver : public IterativeSolver {
  //---------------------------------------------------------------------------------
  // PUBLIC MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
public:
  MultigridSolver(const SparseMatrix   &A,
                  const Vector<double> &b,
                  const unsigned int    pre_nu,
                  const unsigned int    post_nu,
                  const double          tol,
                  const unsigned int    max_iter) :
    IterativeSolver(A, b, tol, max_iter),
    pre_nu(pre_nu), post_nu(post_nu) {
    smoother_tot_iter = 0;
  }

  virtual void
  setup() = 0;

  // Two level multigrid solver
  virtual int
  solve(Vector<double> &x) = 0;

  virtual unsigned int
  get_tot_smoother_iter() const {
    return smoother_tot_iter;
  }

  //---------------------------------------------------------------------------------
  // PROTECTED MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
protected:
  unsigned int pre_nu;
  unsigned int post_nu;

  unsigned int smoother_tot_iter;

  std::unique_ptr<IterativeSolver> pre_smoother;
  std::unique_ptr<IterativeSolver> coarse_smoother;
  std::unique_ptr<IterativeSolver> post_smoother;

  SparseMatrix restrictor;
  SparseMatrix interpolator;

  SparseMatrix A_2h;

  virtual void
  build_restrictor() = 0;

  virtual void
  build_interpolator() = 0;

  virtual void
  coarsen_matrix(SparseMatrix &A) {
    std::cout << "  Coarsening the system matrix..." << std::endl;
    SparseMatrix *A_2h_ = new SparseMatrix;
    *A_2h_              = restrictor.mul(A).mul(interpolator);
    std::cout << "  DONE!" << std::endl;

    A_2h = *A_2h_;
  }
};

#endif