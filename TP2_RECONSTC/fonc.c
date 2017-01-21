#include "fonc.h" 
#include "limace.h"
#include "salade.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>


/* Produit renvoi le produit de => OK !!
 * 2 matrices
 */
Matrix Produit(Matrix Mat1, Matrix Mat2)
{
  int i, j, k;
  double **M1, **M2, **Mres, somme;
  Matrix MatRes;
  
  MatRes = MatCAlloc(Double, MatNbRow(Mat1), MatNbCol(Mat2));
  
  M1 = MatGetDouble(Mat1); 
  M2 = MatGetDouble(Mat2);
  Mres = MatGetDouble(MatRes);
  
  for(i=0; i<MatNbRow(Mat1); i++)
  {
    for(j=0; j<MatNbCol(Mat2); j++)
    {
      somme = 0.0;
      for(k=0; k<MatNbCol(Mat1); k++)
      {
	   somme += M1[i][k] * M2[k][j];
      }
      Mres[i][j] = somme;
    }
  }
  
  return MatRes;
}


/* Transposee renvoi la matrice transposée => OK !!
 * de la matrice mis en paramètre
 */
Matrix Transposee(Matrix Mat)
{
  int i, j;
  double **M, **Mres;
  Matrix MatRes;
  
  MatRes = MatAlloc(Double, MatNbCol(Mat), MatNbRow(Mat));
  
  M = MatGetDouble(Mat); 
  Mres = MatGetDouble(MatRes);
  
  for(i=0; i<MatNbRow(MatRes); i++)
  {
    for(j=0; j<MatNbCol(MatRes); j++)
    {
      Mres[i][j] = M[j][i];
    }
  }
    
  return MatRes;
}

/* renvoi l'indice du minimum 
* d'une matrice M
*/
int min(Matrix M)
{
    int i, j;
    int min_ind;
    double min;
    double **mTab;
    
    mTab = MatGetDouble(M);
    
    min_ind = 0;
    min = mTab[0][0];
    for(i=0; i<MatNbRow(M); i++)
    {
     for(j=0; j<MatNbCol(M); j++)
     {
      if(mTab[i][j] < min)
      {
       min = mTab[i][j];    
       min_ind = i;
      }
     }
    } 
    
    
    return min_ind;
}

int signe(double m34)
{
    int res = (m34<0)?-1:1;    
    return res;
}
