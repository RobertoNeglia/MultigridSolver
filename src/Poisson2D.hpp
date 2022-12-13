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
    virtual double
    value(const double x, const double y) const {};
  };

  class DiffusionCoefficient : public Function2D {
  public:
    double
    value(const double x, const double y) const override {
      return x + y;
    }
  };

  class Gradient {
  public:
    static std::vector<double>
    compute_gradient(const Function2D f, const double x, const double y, const double h) {
      std::vector<double> gradient(DIM);

      std::cout << "=======================" << std::endl;
      std::cout << f.value(x + h, y) << std::endl;

      gradient[0] = (f.value(x + h, y) - f.value(x - h, y)) / (2. * h);
      gradient[1] = (f.value(x, y + h) - f.value(x, y - h)) / (2. * h);

      return gradient;
    }

    static double
    compute_divergence(const Function2D &f, const double x, const double y, const double h) {
      std::vector<double> gradient = compute_gradient(f, x, y, h);
      std::cout << "-----------------------------------" << std::endl;
      std::cout << gradient[0] << gradient[1] << std::endl;
      std::cout << "-----------------------------------" << std::endl;
      return gradient[0] + gradient[1];
    }
  };

  /**
   * Forcing term
   */
  class ForcingTerm : Function2D {
  public:
    double
    value(const double x, const double y) const override {
      return 5;
    }
  };
  /**
   * Dirichlet boundary conditions
   */
  class FunctionG : Function2D {
  public:
    double
    value(const double x, const double y) const override {
      return std::sin(x) + std::cos(y);
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
    A.initialize(m, n);
    generate_discretized_matrix(m, n);

    A.print_matrix();
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

    std::cout << Gradient::compute_divergence(mu, 10, 10, 1.e-6) << std::endl;

      for (unsigned int i = 0; i < m; i++) {
          for (unsigned int j = 0; j < n; j++) {
            A.scalar_mul(Gradient::compute_divergence(mu, i, j, 1.e-6) / (h * h));
          }
      }
  }
};

#endif