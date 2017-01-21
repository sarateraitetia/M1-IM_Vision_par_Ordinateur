/*
*               Fichier contenant toutes les fonctions utiles 
*                       au calibrage d'une cam�ra
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
  int i,j,k;
  
  A = MatCAlloc(Double, 2*MatNbRow(Mat3D), 9); 
  
  Mat3Dtab = MatGetDouble(Mat3D);
  Mat2Dtab = MatGetDouble(Mat2D);
  Atab = MatGetDouble(A);
  
  for(i=0;i < MatNbRow(A);i++)
  {
   if((i%2)==0)
   {           
               j = i/2; /* 1ère ligne */
               k = 0;
   }
   else
   {
              j = (i-1)/2; /* 2ème ligne */
              k = 1;
   }
       
   Atab[i][0+(k*3)] = Mat3Dtab[j][0];
   Atab[i][1+(k*3)] = Mat3Dtab[j][1];	  
   Atab[i][2+(k*3)] = 1;
	
   Atab[i][6] = Mat3Dtab[j][0] * (-Mat2Dtab[j][k]); 
   Atab[i][7] = Mat3Dtab[j][1] * (-Mat2Dtab[j][k]);
   Atab[i][8] = -Mat2Dtab[j][k];
 } 

  return A;
}

            
/* Calibrage 3D : estimation au sens des moindres carres totaux de la matrice
 * de projection perspective à partir de correspondances 3D <--> 2D
 */
Matrix Calibrage(Matrix MatP3, Matrix Matp2)
{    
    
    Matrix A, ata, pVal, pVect, m;
    double **Vect, **M;
    double norme;
    int min_valP;
    int i, j, k;
    
    
    /* Création de la matrice A */
    A = fabA(MatP3, Matp2);
    MatWriteAsc(A, "");
    //afficherMatrix(A);
    
    
    /* Création du vecteur colonne Val (valeur propres) 
    * et de la matrice des vecteurs propres Vect 
    */
    ata = Produit(Transposee(A), A);
    pVal = MatCAlloc(Double, MatNbCol(A), 1); 
    pVect = MatCAlloc(Double, MatNbCol(A),  MatNbCol(A));    
    
    SymEig(ata, &pVal, &pVect);
    //MatWriteAsc(pVect, "");
    //afficherMatrix(pVect);
    
    
    /* Détermination de la valeur propre minimale */
    min_valP = min(pVal);
    printf("min=%d\n\n", min_valP);
    
    
    /* Création du vecteur propre M */
    m = MatCAlloc(Double, 3, 3); 
    M = MatGetDouble(m);
    Vect = MatGetDouble(pVect);
    k=0;
    while(k < MatNbRow(pVect))
    {
         for(i=0; i<MatNbRow(m); i++)
         {
          for(j=0; j<MatNbCol(m); j++)
          {
                           M[i][j] = Vect[k][min_valP];
                           k++;
          }
         }
    }
    
    MatWriteAsc(m, "");

   
   /* Normalisation */
   norme = calcul_norme(M[2], 2);   
   for(i=0; i<MatNbRow(m); i++)
   {
    for(j=0; j<MatNbCol(m); j++)
    {
     M[i][j] = signe(M[MatNbRow(m)-1][MatNbCol(m)-1]) * M[i][j] / norme;
    }
   }
   
   MatWriteAsc(m, "");
   
   
   MatFree(&A);
   MatFree(&ata);
   MatFree(&pVal);
   MatFree(&pVect);
   
   return m;
}
