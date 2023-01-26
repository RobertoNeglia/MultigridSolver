#ifndef _JACOBI_H__
#define _JACOBI_H__

#include "iterative_solver.hpp"

class Jacobi : public IterativeSolver {
  //---------------------------------------------------------------------------------
  // PUBLIC MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
public:
  using IterativeSolver::IterativeSolver;

  int
  solve(Vector<double> &x) override {
    double resid;
    double normb = norm(b);
    if (normb == 0.0)
      normb = 1.0;

    // result of A*x_k
    Vector<double> b_k = A.mul(x);
    // residual at step k: r_k = b - A*x_k = b - b_k;
    Vector<double> r_k(b_k.size());
    r_k = subvec(b, b_k);

      if ((resid = norm(r_k) / normb) <= tol) {
        tol_achieved = resid;
        n_iter       = 0;
        return 0;
    }

    // Preconditioner (Jacobi)
    SparseMatrix   P = get_Jacobi_preconditioner(A);
    Vector<double> p(b_k.size());

      for (unsigned int i = 0; i < max_iter; i++) {
        // calculate the action of the preconditioner on the residual
        p = P.mul(r_k);
        // update the value of the solution
        x = addvec(x, p);
        // update the value of A*x_k
        b_k = A.mul(x);
        // calculate the residual
        r_k = subvec(b, b_k);

          if ((resid = norm(r_k) / normb) <= tol) {
            tol_achieved = resid;
            n_iter       = i;
            return 1;
        }
      }

    tol_achieved = resid;
    n_iter       = max_iter;
    return 2;
  }

  //---------------------------------------------------------------------------------
  // PROTECTED MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
protected:
  SparseMatrix
  get_Jacobi_preconditioner(const SparseMatrix A) const {
    SparseMatrix D;
    D.initialize(A.cols(), A.rows());

      for (unsigned int i = 0; i < A.cols(); i++) {
        D.insert_coeff(1.0 / A.coeff(i, i).first, i, i);
      }

    return D;
  }
};

#endif