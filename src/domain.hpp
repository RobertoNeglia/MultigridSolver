#ifndef _DOMAIN_H__
#define _DOMAIN_H__

#define DIM 2

#include <map>
#include <vector>

#include <Eigen/Eigen>

class Point2D {
private:
  unsigned int x_coord, y_coord;

public:
  Point2D(unsigned int x, unsigned int y) : x_coord(x), y_coord(y) {}

  unsigned int
  scalar_coordinate(unsigned int n_rows) {
    return y_coord * n_rows + x_coord;
  }

  unsigned int
  x() {
    return x_coord;
  }

  unsigned int
  y() {
    return y_coord;
  }
};

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
  unsigned int    n_rows;
  unsigned int    n_cols;
  Eigen::MatrixXd D;

public:
  Domain2D(unsigned int n_rows, unsigned int n_cols) :
    n_rows(n_rows), n_cols(n_cols), D(n_rows, n_cols) {
    // initializes the domain with elements that corresponds to the scalar coordinate inside
    // the matrix
    // for (unsigned int i = 0; i < n_rows; i++) {
    //     for (unsigned int j = 0; j < n_cols; j++) {
    //       D.coeffRef(i, j) = j * n_rows + i;
    //     }
    // }
  }

  std::pair<unsigned int, unsigned int>
  onedim2twodim(unsigned int p) const {
    return std::make_pair(p % n_rows, p / n_rows);
  }

  unsigned int
  twodim2onedim(unsigned int i, unsigned int j) const {
    return j * n_rows + i;
  }

  unsigned int
  rows() const {
    return n_rows;
  }

  unsigned int
  cols() const {
    return n_cols;
  }

  void
  print_domain(std::ostream &os) const {
      for (unsigned int i = 0; i < n_rows; i++) {
          for (unsigned int j = 0; j < n_cols; j++) {
            os << twodim2onedim(i, j) << "\t";
          }
        os << std::endl;
      }
  }
};

#endif