#include <iostream>
#include <vector>

using namespace std;

/*
  TODO
    - inserimento elementi in A in ordine di riga;
    -aggiornare A_col
    -aggiornare A ro3
*/

class Matrix {
private:
  int            nRow;
  int            nCol;
  vector<double> A;
  vector<int>    A_col;
  vector<int>    A_row;

public:
  // commento
  Matrix(int r, int c) : nRow(r), nCol(c) {
      for (unsigned int i = 0; i < nRow; i++) {
        A_row.emplace_back(-1); // -1 indica riga di 0
      }
    A_row.emplace_back(0); //  l'ultimo elemento indica nnz
  }

  void
  insertElement(double val, int i, int j) {
      if (val != 0) {
        A.emplace_back(val);
        A_col.emplace_back(j);
          if (A_row[i] == -1) {
            bool found = false;
              for (unsigned int k = i + 1; k < nRow; k++) {
                  if (!found) {
                      if (A_row[k] != -1) {
                        A_row[i] = A_row[k];
                        found    = true;
                    } else if (A_row[k] != -1)
                      A_row[k]++;
                }
              }
              if (!found) {
                A_row[i] = A_row[A_row.size() - 1];
            }
        }

          for (unsigned int k = i + 1; k < nRow + 1; k++) {
            if (A_row[k] != -1)
              A_row[k]++;
          }
      } else {
          for (unsigned int k = i + 1; k < nRow + 1; k++) {
            if (A_row[k] != -1)
              A_row[k]++;
          }
      }
  }

  double
  getElement(int i, int j) {
    if (A.size() == 0)
      return 0;
    if (A_row[i] == -1)
      return 0;
    bool found = false;
    int  diff;
      for (int k = i + 1; k < nRow + 1 && !found; k++) {
           if (A_row[k] != -1) {
             diff  = A_row[k] - A_row[i];
             found = true;
        }
      }
      if (diff == 1) {
         if (A_col[A_row[i]] == j)
          return A[A_row[i]];
        else
          return 0;
      } else {
           for (int k = 0; k < diff; k++) {
               if ((A_col[A_row[i] + k]) == j) {
                 return A[A_row[i] + k];
            }
          }
        return 0;
      }
  }

  // commento
  void
  print() {
    cout << "Gli elementi sono: " << endl;
      for (unsigned int i = 0; i < A.size(); i++) {
        cout << A[i] << " ";
      }
    cout << endl;
    cout << "Le colonne degli elementi sono: " << endl;
      for (unsigned int i = 0; i < A_col.size(); i++) {
        cout << A_col[i] << " ";
      }
    cout << endl;
    cout << "I primi elementi di ogni riga sono: " << endl;
      for (unsigned int i = 0; i < A_row.size(); i++) {
        cout << A_row[i] << " ";
      }
    cout << endl;

      for (int i = 0; i < nRow; i++) {
          for (int j = 0; j < nCol; j++) {
            std::cout << getElement(i, j) << " ";
          }
        std::cout << std::endl;
      }
    std::cout << std::endl;
  }
};