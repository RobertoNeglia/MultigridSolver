#ifndef _POISSON2D_H__
#define _POISSON2D_H__

#include "Matrix/sparse_matrix_r.hpp"
#include "algebraic_multigrid.hpp"

class Poisson2D {
public:
  class DiffusionCoefficient {
  public:
    double
    value(const double i, const double j) {
      return i + j;
    }
  };

  /**
   * Forcing term
   */
  class ForcingTerm {
  public:
    double
    value(const double i, const double j) {
      return i - j;
    }
  };
  /**
   * Dirichlet boundary conditions
   */
  class FunctionG {
  public:
    double
    value(const double i, const double j) {
      return 1;
    }
  };

  void
  setup() {
    generate_discretized_matrix()
  }

  void
  solve() {
    AlgebraicMultigrid AMG;

    AMG.solve;
  }

private:
  static const int     dim = 2;
  DiffusionCoefficient mu;
  ForcingTerm          f;
  FunctionG            g;
  SparseMatrix         A;

  void
  generate_discretized_matrix(SparseMatrix &A, int m, int n) {
    int mn = m * n;
      for (int i = 0; i < mn; i++) {
        *(A.coeff_ref(i, i).first) = 4;
        if (i > 0 && i % m != 0)
          *(A.coeff_ref(i, i - 1).first) = -1;
        if (i > m - 1)
          *(A.coeff_ref(i, i - m).first) = -1;
        if (i < mn - 1 && (i % m) != (m - 1))
          *(A.coeff_ref(i, i + 1).first) = -1;
        if (i < mn - m)
          *(A.coeff_ref(i, i + m).first) = -1;
      }
  }
};

#endif