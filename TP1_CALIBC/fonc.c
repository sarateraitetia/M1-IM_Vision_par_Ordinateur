/*
*               Fichier contenant toutes les fonctions utiles 
*                       au calibrage d'une caméra
* 
*  par TERAITETIA Sara
*/

#include "fonc.h" 
#include "limace.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

void afficherMatrix(Matrix Mat)
{
  int i, j;
  double **M;
  
  M = MatGetDouble(Mat);
  
  for(i=0; i<MatNbRow(Mat); i++)
  {
    for(j=0; j<MatNbCol(Mat); j++)
    {
      printf("%f\t",M[i][j]);     
      /* printf("M[%d][%d]=%lf\t",i,j,M[i][j]); */
    }
    printf("\n");
  }
}

/* Transposee renvoi la matrice transposée 
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


/* Norme vectorielle 
 *
 */
double calcul_norme(double* vect, int len)
{
    double som=0;
    int i;
    
    for(i=0; i<len; i++)
    {
      som += pow(vect[i],2);
    }
    
    return sqrt(som);
}

/* Fabrication de la matrice A
 * 
 */
Matrix fabA(Matrix Mat3D, Matrix Mat2D)
{
  Matrix A;
  double **Atab, **Mat2Dtab, **Mat3Dtab;
  int i, j;
  
  afficherMatrix(Mat2D);
  afficherMatrix(Mat3D);
  
  A = MatCAlloc(Double, 2*MatNbRow(Mat3D), 12); 
  
  Mat3Dtab = MatGetDouble(Mat3D);
  Mat2Dtab = MatGetDouble(Mat2D);
  Atab = MatGetDouble(A);
  
  for(i=0;i < 2 * MatNbRow(Mat3D);i++){
 
  /* 1ère ligne */
  if((i%2)==0){
	for(j=0;j<3;j++){
	    Atab[i][j] = Mat3Dtab[i/2][j];
     }
	  
	Atab[i][3] = 1;
	
	for(j=8;j<MatNbCol(A)-1;j++){
	  Atab[i][j] = Mat3Dtab[i/2][j-8] * (-Mat2Dtab[i/2][0]); 
	}
	
	Atab[i][MatNbCol(A)-1] = -Mat2Dtab[i/2][0];
 }
 
 /* 2ème ligne */
 else{
	  
   for(j=4;j<7;j++)
   {
    Atab[i][j] = Mat3Dtab[(i-1)/2][j-4]; 
    }
	    
    Atab[i][7] = 1;
	    
    for(j=8;j<MatNbCol(A)-1;j++)
    {
      Atab[i][j] = Mat3Dtab[(i-1)/2][j-8] * (-Mat2Dtab[(i-1)/2][1]); 
    }
	    
    Atab[i][MatNbCol(A)-1] = -Mat2Dtab[(i-1)/2][1];
 }
}
  
  return A;
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
    
    printf("\nVal:\n");
    afficherMatrix(M);
    printf("\n");
    
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

