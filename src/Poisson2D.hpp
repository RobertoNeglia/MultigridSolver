#ifndef _POISSON2D_H__
#define _POISSON2D_H__

#include <cmath>
#include <memory>
#include <vector>

#include "Matrix/sparse_matrix_r.hpp"
#include "algebraic_multigrid.hpp"
#include "domain.hpp"

class Poisson2D {
public:
  Poisson2D(const double x0,
            const double xL,
            const double y0,
            const double yL,
            const double h) :
    x0(x0),
    xL(xL), y0(y0), yL(yL), h(h) {}

  class Function2D {
  public:
    Function2D() {}
    virtual double
    value(const double x, const double y) const = 0;

    virtual double
    value() const = 0;
  };

  class DiffusionCoefficient : public Function2D {
  public:
    DiffusionCoefficient() {}
    double
    value(const double x, const double y) const override {
      return 1.0;
    }

    double
    value() const override {
      return 1.0;
    }
  };

  /**
   * Forcing term
   */
  class ForcingTerm : Function2D {
  public:
    double
    value(const double x, const double y) const override {
      return -5.0 * std::exp(x) * exp(-2.0 * y);
    }

    double
    value() const override {
      return 1.0;
    }
  };
  /**
   * Dirichlet boundary conditions
   */
  class FunctionG : Function2D {
  public:
    double
    value(const double x, const double y) const override {
      return std::exp(x) * std::exp(-2.0 * y);
    }

    double
    value() const override {
      return 1.0;
    }
  };

  class Gradient {
  public:
    static std::vector<double>
    compute_gradient(const Function2D &f, const double x, const double y, const double h) {
      std::vector<double> gradient(DIM);

      gradient[0] = (f.value(x + h, y) - f.value(x - h, y)) / (2. * h);
      gradient[1] = (f.value(x, y + h) - f.value(x, y - h)) / (2. * h);

      return gradient;
    }

    static double
    compute_divergence(const Function2D &f, const double x, const double y, const double h) {
      std::vector<double> gradient = compute_gradient(f, x, y, h);
      return gradient[0] + gradient[1];
    }
  };

  /**
   * Setups all members of the class
   */
  void
  setup() {
    domain.initialize(x0, xL, y0, yL, h);
    domain.print_domain();
    unsigned int m = domain.rows();
    unsigned int n = domain.cols();

    b.resize(m * n);
    generate_rhs(m, n);

    A.initialize(m * n, m * n);
    generate_discretized_matrix(m, n);
    A.print_matrix();

    A.scalar_mul(mu.value());
    A.print_matrix();
  }

  void
  solve() {
    // AlgebraicMultigrid AMG(A, b);

    // AMG.solve(b);
  }

private:
  static const int     dim = 2;
  DiffusionCoefficient mu;
  ForcingTerm          f;
  FunctionG            g;
  SparseMatrix         A;
  std::vector<double>  b;
  double               x0, xL;
  double               y0, yL;
  double               h;
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
    unsigned int mn = m * n;
      for (unsigned int i = 0; i < mn; i++) {
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

    // for (unsigned int i = 0; i < m; i++) {
    //     for (unsigned int j = 0; j < n; j++) {
    //       A.scalar_mul(Gradient::compute_divergence(mu, i, j, 1.e-6) / (h * h));
    //     }
    // }
  }
};

#endif