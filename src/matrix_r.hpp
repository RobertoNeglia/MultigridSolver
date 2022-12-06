#include <iostream>
#include <vector>

class Matrix {
private:
  const int           empty_row = -1;
  int                 nnz;
  int                 n_rows;
  int                 n_cols;
  std::vector<double> A;
  std::vector<int>    A_cols;
  std::vector<int>    A_rows;

  int
  get_end_row(int i) {
      for (int k = i + 1; k < n_rows + 1; k++) {
        if (A_rows[k] != -1)
          return A_rows[k];
      }
    return A_rows[n_rows];
  }

  /*
      0   0   0
      0   0   0
      0   0   0

      A     :
      A_col :
      A_row : -1 -1 -1 0

      insert_in_empty_row(1, 0, 0);
*/

  void
  update_next_A_rows(int i) {
      for (int k = i + 1; i < n_rows; i++) {
        if (A_rows[k] != -1)
          A_rows[k]++;
      }
  }

  void
  update_nnz() {
    A_rows[n_rows]++;
  }

  void
  insert_in_empty_row(double val, int i, int j) {
    int end_row = get_end_row(i);
    A_rows[i]   = end_row;
    update_next_A_rows(i);
    update_nnz();

      if (A.size() == 0) {
        A.emplace_back(val);
        A_cols.emplace_back(j);
    }
  }

public:
  Matrix(int n_rows, int n_cols) : n_rows(n_rows), n_cols(n_cols) {
    nnz = 0;
    for (int i = 0; i < n_rows; i++)
      A_rows.emplace_back(empty_row);
    A_rows.emplace_back(nnz);
  }

  double const
  get_coeff(int i, int j) {
    if (i >= n_rows || j >= n_cols)
      return 0;

    int start_row = A_rows[i];
    if (start_row == -1)
      return 0;

    int end_row = get_end_row(i);
      for (int k = start_row; k < end_row; k++) {
        if (A_cols[k] == j)
          return A[k];
      }
    return 0;
  }

  double &
  get_coeff_ref(int i, int j) {
    int start_row = A_rows[i];
    int end_row   = get_end_row(i);

      for (int k = start_row; k < end_row; k++) {
        if (A_cols[k] == j)
          return A[k];
      }
  }

  void
  insert_coeff(double val, int i, int j) {
      if (val != 0) {
          if (A_rows[i] == -1) { // empty row
            // i have to insert in the exact position in arrays A and A_col, and update
            // A_rows

            insert_in_empty_row(val, i, j);
        }
    }
  }

  void
  print_matrix() {
      for (int i = 0; i < n_rows; i++) {
          for (int j = 0; j < n_cols; j++) {
            std::cout << get_coeff(i, j) << "\t";
          }
        std::cout << std::endl;
      }
  }

  void
  print_structure() {
    std::cout << "A: " << std::endl;
    for (auto i : A)
      std::cout << "\t" << i << " ";
    std::cout << std::endl;

    std::cout << "A_cols: " << std::endl;
    for (auto i : A_cols)
      std::cout << "\t" << i << " ";
    std::cout << std::endl;

    std::cout << "A_rows: " << std::endl;
    for (auto i : A_rows)
      std::cout << "\t" << i << " ";
    std::cout << std::endl;
  }
};