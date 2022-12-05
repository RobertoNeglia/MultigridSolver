#include <vector>
#include <iostream>

using namespace std;


class Matrix {
private:
  int nRow; 
  int nCol;
  vector<double> aa;
  vector<int> ja;
  vector<int> ia; 

public:
  // commento
  Matrix(int r, int c) : nRow(r), nCol(c) {
     
      for(unsigned int i = 0; i < nRow; i++)
         ia.emplace_back(-1);   // -1 indica riga di 0
      ia.emplace_back(0);       //  l'ultimo elemento indica nnz
  }


  void
  insertElement(double val, int i, int j)
  {
    if (val != 0) {
      aa.emplace_back(val);
      ja.emplace_back(j);
      //int howmany = ia[i+1] - ia[i];
      //for(int h = 0; h < howmany; h++)
      //{
        if (ia[i] == -1)
        {
          for (unsigned int k = i+1; k < nRow; k++)
          {
            if (ia[k] != -1)
            {  
              ia[i] = ia[k];
              for ( unsigned j = k; k < nRow; k++)
              {
                ia[k]++;
              }
              break;
            }
            else if (i != 0)
              ia[i] = ia[i-1] + 1;
            else
              ia[i] = 0;  
          }
          
          for (unsigned int k = i+1; k < nRow + 1; k++)
          {
            if (ia[k] != -1)
              ia[k] ++;
          }
        }
        else
        {
          for (unsigned int k = i+1; k < nRow + 1; k++)
          {
            if (ia[k] != -1)
              ia[k] ++;
          }
        }  
      //}
    }

  }
  // commento
  void print() {
    cout << "Gli elementi sono: " << endl;
    for (unsigned int i = 0; i < aa.size(); i++)
    {
      cout << aa[i] << " ";
    }
    cout << endl;
    cout << "Le colonne degli elementi sono: " << endl;
    for (unsigned int i = 0; i < ja.size(); i++)
    {
      cout << ja[i] << " ";
    }
    cout << endl;
    cout << "I primi elementi di ogni riga sono: " << endl;
    for (unsigned int i = 0; i < ia.size(); i++)
    {
      cout << ia[i] << " ";
    }
    cout << endl;
  }
};