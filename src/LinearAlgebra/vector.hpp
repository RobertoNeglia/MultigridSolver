#ifndef _VECTOR_H__
#define _VECTOR_H__

#include <omp.h>
#include <time.h>

#include <cmath>
#include <random>
#include <vector>

template <typename T>
using Vector = std::vector<T>;

template <typename T>
static double
norm(const Vector<T> &a) {
  double n = 0;
#pragma omp parallel for num_threads(omp_get_num_procs())
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
#pragma omp parallel for num_threads(omp_get_num_procs())
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
#pragma omp parallel for num_threads(omp_get_num_procs())
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
#pragma omp parallel for num_threads(omp_get_num_procs())
    for (unsigned int i = 0; i < a.size(); i++) {
      res[i] = a[i] - b[i];
    }
  return;
}

template <typename T>
static void
fill(Vector<T> &v, const T &x) {
  //
#pragma omp parallel for num_threads(omp_get_num_procs())
    for (unsigned int i = 0; i < v.size(); i++) {
      v[i] = x;
    }
}

void
fill_random(Vector<double> &v) {
  std::mt19937                           engine(clock());
  std::uniform_real_distribution<double> dist(-10.0, 10.0);
    for (unsigned int i = 0; i < v.size(); i++) {
      v[i] = dist(engine);
    }
}

void
fill_incremental(Vector<double> &v) {
    for (unsigned int i = 0; i < v.size(); i++) {
      v[i] = i + 1;
    }
}

bool
compare(const Vector<double> &a, const Vector<double> &b, const double &tol) {
  if (a.size() != b.size())
    return false;
  bool eq = true;
    for (unsigned int i = 0; i < a.size(); i++) {
        if (std::abs(a[i] - b[i]) > tol) {
          eq = false;
          break;
      }
    }

  return eq;
}

bool
equal_to(const Vector<double> &v, const double &val) {
  bool eq = true;
    for (unsigned int i = 0; i < v.size(); i++) {
      if (std::abs(val - v[i]) > 1.e-6)
        eq = false;
    }

  return eq;
}

template <typename T>
void
print_vector(const Vector<T> &v) {
    for (T i : v) {
      std::cout << i << " - ";
    }

  std::cout << std::endl;
}

void
press_to_continue(std::ostream &os = std::cout, std::istream &is = std::cin) {
  os << "Press enter to continue: ";
  is.ignore();
}

#endif