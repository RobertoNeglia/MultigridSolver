#ifndef _VECTOR_H__
#define _VECTOR_H__

#include <omp.h>

#include <cmath>
#include <vector>

template <typename T>
using Vector = std::vector<T>;

template <typename T>
static double
norm(const Vector<T> &a) {
  double n = 0;
#pragma omp parallel for num_threads(4)
    for (unsigned int i = 0; i < a.size(); i++) {
      n += a[i] * a[i];
    }
  return sqrt(n);
}

template <typename T>
static Vector<T>
addvec(const Vector<T> &a, const Vector<T> &y) {
  Vector<T> result(a.size(), 0.0);
    if (a.size() != y.size()) {
      std::cout << "VECTORS DO NOT HAVE THE SAME SIZE - ADD" << std::endl;
      return Vector<T>();
  }
#pragma omp parallel for num_threads(4)
    for (unsigned int i = 0; i < a.size(); i++) {
      result[i] = a[i] + y[i];
    }
  return result;
}

template <typename T>
static Vector<T> &
addvec_inplace(Vector<T> &a, const Vector<T> &y) {
    if (a.size() != y.size()) {
      std::cout << "VECTORS DO NOT HAVE THE SAME SIZE - ADD IN PLACE" << std::endl;
      return a;
  }
#pragma omp parallel for num_threads(4)
    for (unsigned int i = 0; i < a.size(); i++) {
      a[i] = a[i] + y[i];
    }
  return a;
}

template <typename T>
static void
subvec(Vector<T> &res, const Vector<T> &a, const Vector<T> &b) {
    if (a.size() != b.size()) {
      std::cout << "VECTORS DO NOT HAVE THE SAME SIZE - SUB" << std::endl;
      return;
  }
#pragma omp parallel for num_threads(4)
    for (unsigned int i = 0; i < a.size(); i++) {
      res[i] = a[i] - b[i];
    }
  return;
}

template <typename T>
static void
            fill(Vector<T> v, const T x) {
#pragma omp parallel for num_threads(4)
    for (unsigned int i = 0; i < v.size(); i++) {
      v[i] = x;
    }
}

template <typename T>
void
print_vector(Vector<T> v) {
    for (T i : v) {
      std::cout << i << " - ";
    }

  std::cout << std::endl;
}

#endif