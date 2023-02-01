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
          if (A.size() == 0) { // empty matrix
            A.emplace_back(val);
            A_col.emplace_back(j);
            A_row[i] = 0;
            A_row[A_row.size() - 1]++;

          } else { // not empty matrix

            // A_row vector update
              if (A_row[i] == -1) { // row i is empty
                bool found = false;
                  for (unsigned int k = i + 1; k < nRow; k++) { // find next row not empty
                      if (!found) {
                          if (A_row[k] != -1) {
                            A_row[i] = A_row[k];
                            found    = true;
                        }
                      } else { // I have found it!!! (update next rows indices)
                        if (A_row[k] != -1)
                          A_row[k]++;
                      }
                  }
                  if (!found) { // update nnz
                    A_row[i] = A_row[A_row.size() - 1];
                }
              } else {
                  for (unsigned int k = i + 1; k < nRow + 1; k++) {
                    if (A_row[k] != -1)
                      A_row[k]++;
                  }
              }
              // update of A, A_col
              if (getElement(i, j) !=
                  0) { // if an element is already written in that position
                getElementRef(i, j) = val;
              } else { // zero in that position
                bool found = false;
                int  diff;                                            // diff è una merda
                  for (int k = i + 1; k < nRow + 1 && !found; k++) { //
                      if (A_row[k] != -1) { // looks for the next non empty row
                        diff =
                          A_row[k] - A_row[i]; // diff is the number of elements on row i
                        found = true;
                    }
                  }
                  if (diff == 1) { // only one element in that line
                    A_col.insert(A_col.begin() + A_row[i], j); // updates A_col (?)
                    A.insert(A.begin() + A_row[i], val);       // updates A (?)

                  } else { // if diff is greater than one (there is the problem that diff
                           // is increased before, while updating A_row)
                    bool found = false;
                      for (int k = 0; k < diff && !found; k++) {
                           if ((A_col[A_row[i] + k]) > j) {
                             A_col.insert(A_col.begin() + A_row[i] + k, j);
                             A.insert(A.begin() + A_row[i] + k, val);
                             found = true;
                        }
                      }
                    if (!found)
                      std::cout << "PIPPO" << endl;
                  }
              }
          }
    }
  }
  // function getElementRef: if present, returns the reference to the element
  double &
  getElementRef(int i, int j) {
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
          return (A[A_row[i]]);
      } else {
           for (int k = 0; k < diff; k++) {
               if ((A_col[A_row[i] + k]) == j) {
                 return A[A_row[i] + k];
            }
          }
      }
  }
  // If present, it returns the element in the given positio
  double
  getElement(int i, int j) {
      if (i < nRow && j < nCol) {
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
    return -1;
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