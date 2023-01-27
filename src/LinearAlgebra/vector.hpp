#ifndef _VECTOR_H__
#define _VECTOR_H__

#include <omp.h>

#include <cmath>
#include <vector>

template <typename T>
using Vector = std::vector<T>;

static double
norm(const Vector<double> &a) {
  double n = 0;
#pragma omp parallel for num_threads(4)
    for (unsigned int i = 0; i < a.size(); i++) {
      n += a[i] * a[i];
    }
  return sqrt(n);
}

static Vector<double>
addvec(const Vector<double> &a, const Vector<double> &y) {
  Vector<double> result(a.size(), 0.0);
    if (a.size() != y.size()) {
      std::cout << "VECTORS DO NOT HAVE THE SAME SIZE" << std::endl;
      return Vector<double>();
  }
#pragma omp parallel for num_threads(4)
    for (unsigned int i = 0; i < a.size(); i++) {
      result[i] = a[i] + y[i];
    }
  return result;
}

static Vector<double> &
addvec_inplace(Vector<double> &a, const Vector<double> &y) {
    if (a.size() != y.size()) {
      std::cout << "VECTORS DO NOT HAVE THE SAME SIZE" << std::endl;
      return a;
  }
#pragma omp parallel for num_threads(4)
    for (unsigned int i = 0; i < a.size(); i++) {
      a[i] = a[i] + y[i];
    }
  return a;
}

static void
subvec(Vector<double> &res, const Vector<double> &a, const Vector<double> &b) {
    if (a.size() != b.size()) {
      std::cout << "VECTORS DO NOT HAVE THE SAME SIZE" << std::endl;
      return;
  }
#pragma omp parallel for num_threads(4)
    for (unsigned int i = 0; i < a.size(); i++) {
      res[i] = a[i] - b[i];
    }
  return;
}

static void
            fill(Vector<double> v, const double x) {
#pragma omp parallel for num_threads(4)
    for (unsigned int i = 0; i < v.size(); i++) {
      v[i] = x;
    }
}

#endif