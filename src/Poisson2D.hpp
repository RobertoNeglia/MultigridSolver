#ifndef _POISSON2D_H__
#define _POISSON2D_H__

#include <vector>

#include "Matrix/sparse_matrix_r.hpp"
#include "algebraic_multigrid.hpp"
#include "domain.hpp"

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
      return i;
    }
  };

  /**
   * Setups all members of the class
   */
  void
  setup() {
    double x0 = 0.0, xL = 1.0;
    double y0 = 0.0, yL = 1.0;
    double h = 1.0 / 4.0;
    domain.initialize(x0, xL, y0, yL, h);
    unsigned int m = domain.rows();
    unsigned int n = domain.cols();
    b.resize(m * n);
    A.initialize(m, n);
    generate_discretized_matrix(m, n);
    generate_rhs(m, n);
  }

  void
  solve() {
    AlgebraicMultigrid AMG(domain, A);

    AMG.solve();
  }

private:
  static const int     dim = 2;
  DiffusionCoefficient mu;
  ForcingTerm          f;
  FunctionG            g;
  SparseMatrix         A;
  std::vector<double>  b;
  Domain2D             domain;

  void
  generate_rhs(const unsigned int m, const unsigned int n) {
      for (unsigned int j = 0; j < n; j++) {
          for (unsigned int i = 0; i < m; i++) {
              if (j == 0 || j == (n - 1) || i == 0 || i == (m - 1)) {
                b[j * m + i] = g.value(i, j);
            } else
              b[j * m + i] = f.value(i, j);
          }
      }
  }

  void
  generate_discretized_matrix(const unsigned int m, const unsigned int n) {
    int mn = m * n;
      for (int i = 0; i < mn; i++) {
        A.insert_coeff(4, i, i);
        if (i > 0 && i % m != 0)
          A.insert_coeff(-1, i, i - 1);
        if (i > m - 1)
          A.insert_coeff(-1, i, i - m);
        if (i < mn - 1 && (i % m) != (m - 1))
          A.insert_coeff(-1, i, i + 1);
        if (i < mn - m)
          A.insert_coeff(-1, i, i + m);
      }
  }
};

#endif