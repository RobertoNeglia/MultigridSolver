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
      ia.emplace_back(0);      //  l'ultimo elemento indica nnz
      
      aa.reserve(nRow);
      ja.reserve(nCol);    
  }


  void
  insertElement(double val, int i, int j)
  { 
    if (val != 0.0) {
      if (ia[i] != 0)
        {
          //aa.emplace_back(val);
          bool alias = false;
          for (unsigned int l = 0; l < NNZRow(i); l++)
          {
            if (ja[ia[i] + l] == j) 
            {
                aa[ia[i] + l] = val;
                alias = true;
                break;
            }    
          }
          if (!alias)
          {
          
            double temp1, temp2;
            int t1, t2;
            temp1 = aa[i];
            t1 = ja[i];
            for (unsigned int n = ia[i] +1 ; n <= aa.size(); n++)
            {   
                temp2 = aa[n];    
                aa[n] = temp1;
                temp1 = temp2;          // inserimento nella posizione corretta
                t2 = ja[n];
                ja[n] = t1;
                t1 = t2;
            } 
            aa[ia[i]] = val;  
            ja[ia[i]] = j;
        //na.emplace_back(n)
        //int howmany = ia[i+1] - ia[i];
        //for(int h = 0; h < howmany; h++)
        //{
        }
        }
        // else if (ia[i] == -1)
        // {
          aa.emplace_back(val);
          ja.emplace_back(j);
          if (i == nRow -1)
            {
              ia[i] = ia[i-1] +1;
              ia[i+1]++;

            }  
          for (unsigned int k = i+1; k < nRow; k++)
          {
            if (ia[k] != -1 )
            {  
              ia[i] = ia[k];
              for ( unsigned n = k; k < nRow; k++)
              {
                ia[k]++;
              }
              break;
            }  
          }
          if (i != 0)
          {
              ia[i] = ia[i-1] + 1;
              ia[nRow]++;
          }
          else
              ia[i] = 0; 
       // }
        // else
        // {
        //   if (NNZRow(i) >= 1)
        //   {
        //     for (unsigned int k = i+1; k < nRow; k++)
        //     {
        //       if (ia[k] != -1)
        //         ia[k] ++;
        //     }
        //     ia[nRow]++;
        //   }
        // }
     }  
      //}
     
    // val = 0: sto inserendo elemento nullo corrisponde ad eliminare un elemento
    else if(aa[i] != 0.0)      
    {
      bool alias = true;
      for (unsigned int l = 0; l < NNZRow(i); l++)
      {
        if (ja[ia[i] + l] == j) 
        {
            aa[ia[i] + l] = val;
            alias = false;
            break;
        }    
      }
      if (alias)
      {
        for (unsigned int n = ia[i] ; n < aa.size(); n++)
        {
            aa[n] = aa[n+1];
            ja[n] = ja[n+1];
        } 
        //aa.erase(aa.size()-1);

      }
      
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


  int NNZRow(int index)
  {
    if (ia[index] == -1)
      return 0;
    for (unsigned int i = 1; i < ia.size() - index; i++)
      if (ia[index + i] != -1)
        return (ia[index+i] - ia[index]);
    return 0;
  }

  // void sovrascrivi(double val1, double val2)
  // {

  // }

  // bool binarySearch(double valore, vector<double> a, int inizio, int fine)
  // {
      
  // }
};

