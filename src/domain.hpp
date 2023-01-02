#ifndef _DOMAIN_H__
#define _DOMAIN_H__

#define DIM 2

/**
 * Domain is composed of elements u_i, that are the solution to the linear system that
 * will be constructed A * u = f.
 *
 * Matrix A comes from the discretization of the 2D poisson problem using finite
 * differences and a five-point Laplacian stencil.
 *
 * Vector f is the forcing term of our problem, the independent known data
 *
 * Vector u is the solution of the linear system that is the evaluation
 */
class Domain2D {
private:
  double       x0, xL;
  double       y0, yL;
  double       h;
  unsigned int n, m;

public:
  Domain2D() {}

  void
  initialize(double x0_, double xL_, double y0_, double yL_, double h_) {
    x0 = x0_, xL = xL_, y0 = y0_, yL = yL_, h = h_;
    n = (xL - x0) / h;
    m = (yL - y0) / h;
  }

  static std::pair<unsigned int, unsigned int>
  onedim2twodim(unsigned int p, unsigned int _m) {
    return std::make_pair(p % _m, p / _m);
  }

  static unsigned int
  twodim2onedim(unsigned int i, unsigned int j, unsigned int m) {
    return j * m + i;
  }

  unsigned int
  rows() const {
    return m;
  }

  unsigned int
  cols() const {
    return n;
  }

  void
  print_domain(std::ostream &os = std::cout) const {
      for (unsigned int i = 0; i < m; i++) {
          for (unsigned int j = 0; j < n; j++) {
            os << twodim2onedim(i, j, m) << "\t";
          }
        os << std::endl;
      }
  }
};

#endif