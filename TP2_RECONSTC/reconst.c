#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "limace.h"
#include "salade.h"
#include "fonc.h"


/*----------------------------------------------------------------------
* ----------------------- Programme Principal --------------------------
*-----------------------------------------------------------------------
*/
int main(int argc, char *argv[])
{
    int i,j,k;
    Matrix p2d_gauche, p2d_droite, Mg, Md, MatA;
	Matrix ata, pVal, pVect, m;
	int min_valP;
    Matrix p3d_reconst=NULL;
    double **A, **MatG, **MatD, **MatPG, **MatPD; 
	double **Vect, **M;
	double **p3d;
    
    if(argc!=6)
    {
     printf("Usage: ./reconst <Mg> <Md> <p2d_gauche> <p2d_droite> <p3d_reconst>\n");
     return 1;
    }
    
	Mg = MatReadAsc(argv[1]);
	Md = MatReadAsc(argv[2]);
    p2d_gauche = MatReadAsc(argv[3]);
    p2d_droite = MatReadAsc(argv[4]);
	
    MatA = MatCAlloc(Double, 4, 4); 
	
    A = MatGetDouble(MatA);
	MatG = MatGetDouble(Mg);
	MatD = MatGetDouble(Md);
	MatPG = MatGetDouble(p2d_gauche); 
	MatPD = MatGetDouble(p2d_droite);
  
    // Construction de la matrice p3d_reconst
	p3d_reconst = MatCAlloc(Double, 80, 3);
	p3d = MatGetDouble(p3d_reconst);
	for(j=0; j<80; j++)
	{
		for(i=0;i<3;i++)
		{
			A[0][i] = MatG[0][i] - MatPG[j][0] * MatG[2][i];
			A[1][i] = MatG[1][i] - MatPG[j][1] * MatG[2][i];
			A[2][i] = MatD[0][i] - MatPD[j][0] * MatD[2][i];
			A[3][i] = MatD[1][i] - MatPD[j][1] * MatD[2][i];
		}
    
		A[0][i] = -(MatPG[j][0] * MatG[2][i] - MatG[0][i]);
		A[1][i] = -(MatPG[j][1] * MatG[2][i] - MatG[1][i]);
		A[2][i] = -(MatPD[j][0] * MatD[2][i] - MatD[0][i]);
		A[3][i] = -(MatPD[j][1] * MatD[2][i] - MatD[1][i]);		
		
		// Ce qui suit :
		// - SymEigh de A
		// - on prend le vecteur propre associé à la plus petite valeur propre
		// - on normalise : on divise chaque valeur par la dernière valeur du vecteur B
		
		// Création du vecteur colonne Val (valeur propres) 
		// et de la matrice des vecteurs propres Vect 
		ata = Produit(Transposee(MatA), MatA);
		pVal = MatCAlloc(Double, MatNbCol(MatA), 1); 
		pVect = MatCAlloc(Double, MatNbCol(MatA),  MatNbCol(MatA));    
    
		SymEig(ata, &pVal, &pVect);
    
		// Détermination de la valeur propre minimale 
		min_valP = min(pVal);    
    
		// Création du vecteur propre M 
		m = MatCAlloc(Double, 4, 1); 
		M = MatGetDouble(m);
		Vect = MatGetDouble(pVect);
		k=0;
		for(i=0; i<MatNbRow(m); i++)
		{
			M[i][0] = Vect[k][min_valP];
            k++;
		}
		
		// Normalisation et Remplissage de la matrice "p3d_reconst"
		for(i=0; i<MatNbCol(p3d_reconst); i++)
		{
			p3d[j][i] = M[i][0] / M[3][0];
		}
	}
	
	// Ecriture de la matrice reconstruite dans le fichier
	MatWriteAsc(p3d_reconst,argv[5]);
	
    
	/* Liberation des matrices */
	MatFree(&p2d_gauche); 
    MatFree(&p2d_droite);
    MatFree(&Mg);
    MatFree(&Md);
    MatFree(&MatA);
    MatFree(&ata);
    MatFree(&pVal);
    MatFree(&pVect);
    MatFree(&m); 
	MatFree(&p3d_reconst); 
	
    return 0;
}










