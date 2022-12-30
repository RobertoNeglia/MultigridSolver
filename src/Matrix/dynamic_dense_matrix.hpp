#ifndef _DYNAMIC_DENSE_MATRIX_H__
#define _DYNAMIC_DENSE_MATRIX_H__

#include <iostream>
#include <map>
#include <vector>

class MatrixDI {
  //---------------------------------------------------------------------------------
  // PRIVATE MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
private:
  unsigned int                    zero = 0;
  unsigned int                    m = 0, n = 0;
  std::vector<std::map<int, int>> M;

  //---------------------------------------------------------------------------------
  // PUBLIC MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
public:
  MatrixDI() {
    zero = 0;
    m    = 0;
    n    = 0;
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
  insert_coeff(const int val, const unsigned int i, const unsigned int j) {
      if (i + 1 > M.size()) {
        M.resize(i + 1);
        m = i + 1;
    }

    const auto it = M[i].find(j);
      if (it == M[i].end()) { // element not present
        n                                     = std::max(n, j + 1);
        (*M.at(i).emplace(j, 0).first).second = val;
    }
    (*it).second = val;
  }

  // element access in read only (const version)
  const int &
  coeff(const unsigned int i, const unsigned int j) const {
    return M.at(i).at(j);
  }

  // element access in write (returns non-const reference)
  int &
  coeff_ref(const unsigned int i, const unsigned int j) {
    return M.at(i).at(j);
  }

  void
  print(std::ostream &os = std::cout) {
      for (unsigned int i = 0; i < M.size(); i++) {
          for (auto [j, key] : M.at(i)) {
            os << "(" << i << "," << j << ") = " << key << std::endl;
          }
      }
  }
};

#endif