#ifndef _JACOBI_H__
#define _JACOBI_H__

#include "iterative_solver.hpp"

class Jacobi : public IterativeSolver {
  //---------------------------------------------------------------------------------
  // PUBLIC MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
public:
  Jacobi(const SparseMatrix   &A,
         const Vector<double> &b,
         const double          tol,
         const unsigned int    max_iter) :
    IterativeSolver(A, b, tol, max_iter) {
    P = get_Jacobi_preconditioner(A);
  }

  int
  solve(Vector<double> &x) override {
    double resid;
    double normb = norm(b);
    if (normb == 0.0)
      normb = 1.0;

    // result of A*x_k
    Vector<double> b_k(x.size());
    A.mul(b_k, x);

    // residual at step k: r_k = b - A*x_k = b - b_k;
    Vector<double> r_k(b_k.size());
    subvec(r_k, b, b_k);

      if ((resid = norm(r_k) / normb) <= tol) {
        tol_achieved = resid;
        n_iter       = 0;
        return 0;
    }

    Vector<double> p(b_k.size());

      for (unsigned int i = 0; i < max_iter; i++) {
        // calculate the action of the preconditioner on the residual
        P.mul(p, r_k);

        // update the value of the solution
        addvec_inplace(x, p);

        // update the value of A*x_k
        A.mul(b_k, x);

        // calculate the residual
        subvec(r_k, b, b_k);

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
  SparseMatrix P;

  SparseMatrix &
  get_Jacobi_preconditioner(const SparseMatrix A) const {
    SparseMatrix *D = new SparseMatrix;
    D->initialize(A.cols(), A.rows());

      for (unsigned int i = 0; i < A.cols(); i++) {
        D->insert_coeff(1.0 / A.coeff(i, i).first, i, i);
      }

    return *D;
  }
};

#endif