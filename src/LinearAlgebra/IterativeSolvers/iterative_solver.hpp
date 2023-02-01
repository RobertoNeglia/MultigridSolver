#ifndef _ITERATIVE_SOLVER_H__
#define _ITERATIVE_SOLVER_H__

#include "../sparse_matrix.hpp"

// Abstract iterative solver class
class IterativeSolver {
  //---------------------------------------------------------------------------------
  // PUBLIC MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
public:
  IterativeSolver(const SparseMatrix   &A,
                  const Vector<double> &b,
                  const double         &tol,
                  const unsigned int   &max_iter) :
    A(A),
    b(b), tol(tol), max_iter(max_iter) {
    tol_achieved = 0;
    n_iter       = 0;
  }

  virtual int
  solve(Vector<double> &x) = 0;

  void
  change_rhs(const Vector<double> &new_b) {
    b = new_b;
  }

  unsigned int
  get_iter() const {
    return n_iter;
  }

  double
  get_tol_achieved() const {
    return tol_achieved;
  }

  //---------------------------------------------------------------------------------
  // PROTECTED MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
protected:
  SparseMatrix   A;
  Vector<double> b;

  const double       tol;
  const unsigned int max_iter;
  double             tol_achieved;
  unsigned int       n_iter;
};

#endif