#include <algorithm>
#include <cmath>
#include <vector>

#include "matrix.hpp"

using namespace std;

class Jacobi {
private:
  Matrix              A;
  std::vector<double> b;
  double              tol;
  int                 max_iter;

public:
  Jacobi(Matrix A, std::vector<double> b) : A(A), b(b) {
    tol      = 1.e-6;
    max_iter = 1000;
  }

  Jacobi(Matrix A, std::vector<double> b, double tol, int max_iter) :
    A(A), b(b), tol(tol), max_iter(max_iter) {}

  void
  solve(std::vector<double> &x) {
    Matrix P = Matrix(A.getRows(), A.getCols());
    // cout << P.getCols() << " " << P.getRows() << endl;
    std::vector<double> a = A.Multiplication(x);

    std::vector<double> r(a.size(), 0.0);

    // cout << r.size() << endl;

    P = A.JacPrec();
    // cout << r.size() << endl;
    std::vector<double> p(a.size());

      for (unsigned int i = 0; i < max_iter; i++) {
        r = minvec(b, a);
        p = P.Multiplication(r);
          // for (unsigned int i = 0; i < a.size(); i++)
          //   {
          //     cout << r[i] << " ";
          //   }
          //   cout << endl;

          // cout << norm(r) / norm(b) << endl;
          if (norm(r) / norm(b) <= tol) {
            x = addvec(x, p);
            break;
        }

        // cout << "lele";
        x = addvec(x, p);
        // for (unsigned int i = 0; i < a.size(); i++)
        //   {
        //     cout << x[i] << " ";
        //   }

        //   //cout << endl;
        // cout << "Sono A * x" << endl;
        a = A.Multiplication(x);
        // for (unsigned int i = 0; i < a.size(); i++)
        //   {
        //     cout << r[i] << " ";
        //   }
        // cout << "Sono P * r" << endl;
        // //p = P.Multiplication(r);
        // if (i < 5)
        // {
        // for (unsigned int j = 0; j < a.size(); j++)
        //   {
        //     cout << a[j] << " ";
        //   }
        // }
      }
    // return x;
  }

  void
  setTol(double tol) {
    tol = tol;
  }
  double
  norm(std::vector<double> a) {
    double n = 0;
      for (unsigned int i = 0; i < a.size(); i++) {
        n += a[i] * a[i];
      }
    return sqrt(n);
  }

  vector<double>
  addvec(vector<double> a, vector<double> y) {
    vector<double> result(a.size(), 0.0);
      if (a.size() != y.size()) {
        cout << "lele non e' possibile";
        exit(0);
    }
      // result.reserve(a.size());
      for (unsigned int i = 0; i < a.size(); i++) {
        result[i] = a[i] + y[i];
      }
    return result;
  }

  vector<double>
  minvec(vector<double> a, vector<double> y) {
    vector<double> result(a.size(), 0.0);
      if (a.size() != y.size()) {
        cout << "lele non e' possibile";
        exit(0);
    }

      // result.reserve(a.size());
      for (unsigned int i = 0; i < a.size(); i++) {
        result[i] = a[i] - y[i];
      }
    // cout << result.size() << endl;
    return result;
  }
};
