#ifndef __SPARSE_MATRIX_R_H__
#define __SPARSE_MATRIX_R_H__

#include <iostream>
#include <map>
#include <vector>

/*
  TODO:
    - operator access to elements of the matrix
*/

enum Increment { DECREASE = -1, INCREASE = 1 };

/**
 * Implementation of sparse matrix that uses the CSR (Compressed Sparse Row) format to
 * store the sparse data
 */
class SparseMatrix {
  //---------------------------------------------------------------------------------
  // PRIVATE MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
private:
  // Matrix parameters
  const int    empty_row = -1;
  int          nnz;
  unsigned int n_rows;
  unsigned int n_cols;

  // Storing data structures
  std::vector<double>       A;
  std::vector<unsigned int> A_col;
  std::vector<int>          A_row;

  /**
   * Returns the end of row i of the matrix inside the vectors A and A_col
   */
  int
  end_row(const unsigned int i) const {
      for (unsigned int k = i + 1; k < n_rows + 1; k++) {
        if (A_row[k] != -1)
          return A_row[k];
      }
    return A_row[n_rows];
  }

  /**
   * Updates all elements after the current in the vector A_row.
   *
   * If an element is inserted, Increment should be INCREASE, otherwise if an element is
   * removed (= 0 is inserted in that position), DECREASE should be used
   */
  void
  update_next_A_row(const unsigned int i, const Increment incr) {
      for (unsigned int k = i + 1; k < n_rows; k++) {
          if (A_row[k] != -1) {
            A_row[k] += incr;
        }
      }
  }

  /**
   * Updates the number of nnz elements inside the matrix, both in the A_row vector and in
   * Matrix class variable
   */
  void
  update_nnz(const Increment incr) {
    nnz += incr;
    A_row[n_rows] += incr;
  }

  /**
   * Finds the position inside vectors A and A_cols of the element with coordinates (i,j)
   */
  unsigned int
  find_position(const unsigned int i, const unsigned int j) const {
    const unsigned int start = A_row[i];
    const unsigned int end   = end_row(i);

    unsigned int k;
    unsigned int col;
      for (k = start; k < end; k++) {
        col = A_col[k];
        if (j < col)
          return k;
      }
    return k;
  }

  /**
   * Inserts an elment in an empty row
   */
  void
  insert_in_empty_row(const double val, const int i, const int j) {
    int end = end_row(i);

    // A and A_col vectors update
      if (A.size() == 0) { // structures are empty
        A.emplace_back(val);
        A_col.emplace_back(j);
      } else {
        A.insert(A.begin() + end, val);
        A_col.insert(A_col.begin() + end, j);
      }
    // end of A and A_col vectors update

    // A_row vector update
    A_row[i] = end;
    update_next_A_row(i, INCREASE);
    update_nnz(INCREASE);
    // end of A_row vector update
  }

  /**
   * Inserts an element in a non-empty row
   */
  void
  insert_in_non_empty_row(const double val, const int i, const int j) {
    std::pair<double, bool> el = coeff(i, j);
      if (el.second && el.first == 0) { // if element is not already present in the matrix

        // A and A_col vectors update
        int k = find_position(i, j);
        A.insert(A.begin() + k, val);
        A_col.insert(A_col.begin() + k, j);
        // end of A and A_col vectors update

        // A_row vector update
        update_next_A_row(i, INCREASE);
        update_nnz(INCREASE);
        // end of A_row vector update

      } else if (el.second && el.first != 0) {
        // if element is present in the matrix, simply override it
        std::pair<double *, bool> cf = coeff_ref(i, j);
        if (cf.second)
          *(cf.first) = val;
    }
  }

  /**
   * Removes an element (= 0 is inserted in that position) from the matrix
   */
  void
  remove(const unsigned int i, const unsigned int j) {
    std::pair<double, bool> el = coeff(i, j);
      // check that the element is present in the matrix
      if (el.second && el.first != 0) { // element is in the matrix
        const unsigned int start = A_row[i];
        const unsigned int end   = end_row(i);

        // A and A_col vectors update
        int k = find_position(i, j) - 1;
        A.erase(A.begin() + k);
        A_col.erase(A_col.begin() + k);
        // end of A and A_col vectors update

        // A_row vector update
        const int diff = end - start;
        if (diff == 1)
          A_row[i] = -1;
        update_next_A_row(i, DECREASE);
        update_nnz(DECREASE);
        // end of A_row vector update
    }
  }

  //---------------------------------------------------------------------------------
  // PUBLIC MEMBERS DECLARATION
  //---------------------------------------------------------------------------------
public:
  /**
   * Matrix class constructor
   */
  SparseMatrix() {}

  void
  initialize(const int m, const int n) {
    // Initialize variables and structures
    n_rows = m;
    n_cols = n;
    nnz    = 0;
    for (unsigned int i = 0; i < n_rows; i++)
      A_row.emplace_back(empty_row);
    A_row.emplace_back(nnz);
  }

  unsigned int
  rows() const {
    return n_rows;
  }

  unsigned int
  cols() const {
    return n_cols;
  }

  unsigned int
  get_nnz() const {
    return nnz;
  }

  /**
   * Returns the reference to the element, if present
   */
  std::pair<double *, bool>
  coeff_ref(const unsigned int i, const unsigned int j) {
    int start = A_row[i];   // beginning of the row i in the vectors A and A_col
    int end   = end_row(i); // end of the row i in the vectors A and A_col

      for (int k = start; k < end; k++) {
        if (A_col[k] == j) // if the element is present, return it
          return std::make_pair(&A[k], true);
      }
    return std::make_pair(nullptr, false);
  }

  /**
   * Returns the coefficient in position (i,j) of the matrix
   */
  std::pair<double, bool> const
  coeff(const unsigned int i, const unsigned int j) const {
    if (i >= n_rows || j >= n_cols) // out of bound requests, return false
      return std::make_pair(0, false);

    int start = A_row[i];
    if (start == -1) // empty row, returns 0
      return std::make_pair(0, true);

    // not empty row
    int end = end_row(i); // get end of the row
      for (int k = start; k < end; k++) {
        if (A_col[k] == j) // if element is found, return it
          return std::make_pair(A[k], true);
      }
    return std::make_pair(0, true); // otherwise always returns 0
  }

  /**
   * Inserts an elements inside the matrix at position (i,j)
   */
  void
  insert_coeff(const double val, const int i, const int j) {
      if (val != 0) {
        // I have to insert in the exact position in arrays A and A_col, and update
        // A_row
          if (A_row[i] == -1) { // empty row
            insert_in_empty_row(val, i, j);
          } else { // not empty row
            insert_in_non_empty_row(val, i, j);
          }
      } else {
        // I have to remove an element
        remove(i, j);
      }
  }

  // computes A*v
  std::pair<std::vector<double>, bool>
  mul(const std::vector<double> &v) const {
    if (n_cols != v.size())
      return std::make_pair(std::vector<double>(), false);

    std::vector<double> Av(v.size());

      for (unsigned int i = 0; i < n_rows; i++) {
        double sum   = 0.0;
        int    start = A_row[i];
        int    end   = end_row(i);

        for (int k = start; k < end; k++)
          sum += A[k] * v[A_col[k]];

        Av[i] = sum;
      }

    return std::make_pair(Av, true);
  }

  void
  scalar_mul(const double alpha) {
      for (unsigned int i = 0; i < A.size(); i++) {
        A[i] *= alpha;
      }
  }

  /**
   * Prints the matrix in a human readable way
   */
  void
  print_matrix(std::ostream &os = std::cout) const {
    os << "Number of nnz elements: " << nnz << std::endl;
      for (unsigned int i = 0; i < n_rows; i++) {
          for (unsigned int j = 0; j < n_cols; j++) {
            os << coeff(i, j).first << "\t";
          }
        os << std::endl;
      }
  }

  /**
   * Prints the matrix as it is stored
   */
  void
  print_structure(std::ostream &os = std::cout) const {
    os << "A: " << std::endl;
    for (auto i : A)
      os << "\t" << i << " ";
    os << std::endl;

    os << "A_col: " << std::endl;
    for (auto i : A_col)
      os << "\t" << i << " ";
    os << std::endl;

    os << "A_row: " << std::endl;
    for (auto i : A_row)
      os << "\t" << i << " ";
    os << std::endl;
  }
};

#endif