/*
    Questa struttura funge da efficient matrix container.
    AA  Array degli elementi non zero
    JA  Array degli indici di colonna degli elementi non zero
    IA  Array col numero di nnz elements fino all'i-esima riga

    Si inserisce con insertionElement
    Per eliminare un elemento inserire 0 in quella posizione.
    Se si inserisce un valore in una cella piena, quest'ultimo sostituisce il precedente.
    Se si inserisce 0 in una cella vuota non si ha alcun effetto.

    TODO: funzione MatrixMul
          Jacobi
*/

#include <iostream>
#include <vector>

using namespace std;

class Matrix {
private:
  int            nRow;
  int            nCol;
  vector<double> aa;
  vector<int>    ja;
  vector<int>    ia;

public:
  // commento
  Matrix(int r, int c) : nRow(r), nCol(c) {
    for (unsigned int i = 0; i < nRow; i++)
      ia.emplace_back(-1); // -1 indica riga di 0
    ia.emplace_back(0);    //  l'ultimo elemento indica nnz

    aa.reserve(nRow);
    ja.reserve(nCol);
  }

  int
  getRows() {
    return nRow;
  }

  int
  getCols() {
    return nCol;
  }

  Matrix
  JacPrec() {
    Matrix w(nRow, nCol);
      for (unsigned int j = 0; j < nRow; j++) {
          for (unsigned int i = ia[j]; i < ia[j] + NNZRow(j); i++) {
              if (j == ja[i]) {
                w.aa.emplace_back(1 / aa[i]);
                w.ja.emplace_back(j);
                w.ia[j] = w.ia[nRow];
                w.ia[nRow]++;
                break;
            }
          }
      }
    return w;
  }

  void
  insertElement(double val, int i, int j) {
      if (val != 0.0) {
          if (ia[i] != -1) {
            // aa.emplace_back(val);
            bool alias = false;
              for (unsigned int l = 0; l < NNZRow(i); l++) {
                  if (ja[ia[i] + l] == j) {
                    aa[ia[i] + l] = val;
                    alias         = true;
                    break;
                }
              }
              if (!alias) {
                int    start = 0;
                double temp1, temp2;
                int    t1, t2;

                aa.emplace_back(-4);
                ja.emplace_back(-4);

                  for (unsigned int m = ia[i]; m < ia[i] + NNZRow(i); m++) {
                      if (j < ja[m]) {
                        start = m;
                        temp1 = aa[m];
                        t1    = ja[m];
                          // int j = 0;
                          // while (j > start)
                          // {
                          //   start = ja[ja[i] + j];
                          //   j++;
                          // }
                          for (unsigned int n = start + 1; n < aa.size(); n++) {
                            temp2 = aa[n];
                            aa[n] = temp1;
                            temp1 = temp2; // inserimento nella posizione corretta
                            t2    = ja[n];
                            ja[n] = t1;
                            t1    = t2;
                          }
                        break;
                      } else if (m == ia[i] + NNZRow(i) - 1) {
                        start = m + 1;
                        temp1 = aa[m + 1];
                        t1    = ja[m + 1];
                          for (unsigned int n = start + 1; n < aa.size(); n++) {
                            temp2 = aa[n];
                            aa[n] = temp1;
                            temp1 = temp2; // inserimento nella posizione corretta
                            t2    = ja[n];
                            ja[n] = t1;
                            t1    = t2;
                          }
                        break;
                    }
                  }
                aa[start] = val;
                ja[start] = j;
                  for (unsigned int k = i + 1; k < nRow; k++) {
                    if (ia[k] != -1)
                      ia[k]++;
                  }
                ia[nRow]++;
                // na.emplace_back(n)
                // int howmany = ia[i+1] - ia[i];
                // for(int h = 0; h < howmany; h++)
                //{
            }
          } else if (ia[i] == -1) {
            int    start   = 0;
            int    aaindex = 0;
            double p1, p2;
            int    t1, t2;

              for (unsigned int v = i + 1; v < ia.size(); v++) {
                  if (ia[v] != -1) {
                    start   = v;
                    aaindex = ia[v];
                    break;
                }
              }
            p1 = aa[aaindex];
            t1 = ja[aaindex];
            aa.emplace_back(-4);
            ja.emplace_back(-4);
              for (unsigned int n = aaindex + 1; n < aa.size(); n++) {
                p2    = aa[n];
                aa[n] = p1;
                p1    = p2; // inserimento nella posizione corretta
                t2    = ja[n];
                ja[n] = t1;
                t1    = t2;
              }

            aa[aaindex] = val;
            ja[aaindex] = j;
              // aa.emplace_back(val);
              // ja.emplace_back(j);
              if (i == nRow - 1) {
                ia[i] = ia[nRow];
                ia[nRow]++;
              } else if (i != 0 && i != nRow - 1) {
                // ia[i] = ia[i-1] + 1;
                ia[i] = ia[start];
                  for (unsigned int k = start; k < nRow; k++) {
                      if (ia[k] != -1) {
                        ia[k]++;
                    }
                  }
                ia[nRow]++;
              } else if (i == 0) {
                ia[i] = 0;
                ia[nRow]++;
            }
        }

      }
      //}

      // val = 0: sto inserendo elemento nullo corrisponde ad eliminare un elemento
      else {
        bool alias = false;
          for (unsigned int l = 0; l < NNZRow(i); l++) {
              if (ja[ia[i] + l] == j) {
                // aa[ia[i] + l] = val;
                alias = true;
                break;
            }
            // alias = false;
          }
          if (alias) {
              for (unsigned int n = ia[i]; n < aa.size(); n++) {
                aa[n] = aa[n + 1];
                ja[n] = ja[n + 1];
              }
            // aa.erase(aa.size()-1);
            aa.pop_back();
            ja.pop_back();
              for (unsigned int k = nRow; k > i; k--) {
                if (ia[k] != -1)
                  ia[k]--;
              }
            if (NNZRow(i) == 0)
              ia[i] = -1;
        }
      }
  }
  // commento
  void
  print() {
    cout << "Gli elementi sono: " << endl;
      for (unsigned int i = 0; i < aa.size(); i++) {
        cout << aa[i] << " ";
      }
    cout << endl;
    cout << "Le colonne degli elementi sono: " << endl;
      for (unsigned int i = 0; i < ja.size(); i++) {
        cout << ja[i] << " ";
      }
    cout << endl;
    cout << "I primi elementi di ogni riga sono: " << endl;
      for (unsigned int i = 0; i < ia.size(); i++) {
        cout << ia[i] << " ";
      }
    cout << endl;
  }

  int
  NNZRow(int index) {
    int nnz = 0;
    if (ia[index] == -1)
      return nnz;
      for (unsigned int i = index + 1; i < ia.size(); i++) {
          if (ia[i] != -1) {
            nnz = ia[i] - ia[index];
            break;
        }
      }
    return nnz;
    // if (ia[index] == -1)
    //   return 0;
    // for (unsigned int i = 1; i < ia.size() - index; i++)
    //   if (ia[index + i] != -1)
    //     return (ia[index+i] - ia[index]);
    // return 0;
  }

  vector<double>
  Multiplication(vector<double> x) {
    vector<double> b(nRow, 0.0);
    // b.reserve(nRow);
      if (nCol != x.size()) {
        cout << "Problem not well-posed" << endl;
        exit(0);
    }
      for (unsigned int i = 0; i < nRow; i++) {
          // b[i] = 0;
          for (unsigned int j = ia[i]; j < ia[i] + NNZRow(i); j++) {
            b[i] += aa[j] * x[ja[j]];
            // cout << "ITERAZIONE NUMERO: " << i << endl << endl;
            // cout << "Giro numero: " << j - ia[i];
            // cout << "Valore: " << b[i] << endl;
          }
        // for (unsigned int i = 0; i < nRow; i++)
        //   {
        //     cout << b[i] << " ";
        //   }
      }
    return b;
  }
  // void sovrascrivi(double val1, double val2)
  // {

  // }

  // bool binarySearch(double valore, vector<double> a, int inizio, int fine)
  // {

  // }
};
