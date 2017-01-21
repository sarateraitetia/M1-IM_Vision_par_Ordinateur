#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "fonc.h"
#include "erreurs.h"

/* Codes de retour */
#define OK 0
#define HELP 1
#define ERR_NB_PARAM 2
#define ERR_MAT_3D 3
#define ERR_MAT_3D_TYPE 4
#define ERR_MAT_2D 5
#define ERR_MAT_2D_TYPE 6
#define ERR_MAT_OUT 7


void Syntaxe(void)
{
  Usage("MatricePoints3D MatricePoints2D MatriceProjectionPerspective\n"
        "-h\n");
}

void Description(void)
{
  Mesg("ROLE\n");
  Mesg("\tCalibrage 3D d'une camera\n");
  Mesg("ARGUMENTS\n");
  Mesg("\tMatricePoints3D : Matrice nx3 (format Matrix) des coordonnees des points de la mire\n");
  Mesg("\tMatricePoints2D : Matrice nx2 (format Matrix) des coordonnees des projections dans l'image des points de la mire\n");
  Mesg("\tMatriceProjectionPerspective : nom du fichier destination qui va contenir la matrice 3x4 (format Matrix) de projection perspective\n");
  Mesg("OPTION\n");
  Mesg("\t-h : affichage de l'aide\n");
  Mesg("DIAGNOSTIC (codes de retour)\n");
  Mesg("\t0 : operation realisee sans probleme\n");
  Mesg("\t1 : aide demandee\n");
  Mesg("\t2 : mauvais nombre de parametres\n");
  Mesg("\t3 : probleme d'ouverture du fichier des coordonnees 3D\n");
  Mesg("\t4 : type de matrice 3D incorrect (Double attendu)\n");
  Mesg("\t5 : probleme d'ouverture du fichier des coordonnees 2D\n");
  Mesg("\t6 : type de matrice 2D incorrect (Double attendu)\n");
  Mesg("\t7 : probleme lors du calcul de la matrice de projection perspective\n");
}

int main(int argc, char *argv[])
{
  /* Declaration de la fonction qui sera definie apres la fonction main */
  Matrix Calibrage(Matrix MatP3, Matrix Matp2);
  /* Declaration des variables de la fonction main */
  Matrix MatP3=NULL,Matp2=NULL,MatProjPers=NULL;
  
  /* Initialisation du mecanisme d'affichage des messages */
  InitMesg(argv);
  /* Verification du nombre de parametres */
	if (argc!=2 && argc!=4)
	{
	  Syntaxe();
	  return ERR_NB_PARAM;
	}
	if (argc==2)
	{
	  if (!strcmp(argv[1],"-h"))
	  {
	    Syntaxe();
	    Description();
	    return HELP;
	  }
	  else /* un seul parametre different de -h */
	  {
	    Syntaxe();
	    return ERR_NB_PARAM;
	  }
	}
	/* Lecture de la matrice des points 3D */
	MatP3=MatReadAsc(argv[1]);
	if (MatP3==NULL)
	{
		Erreur("Probleme lors de la lecture de la matrice des coordonnees 3D");
		return ERR_MAT_3D;
	}
	/* Verification du type de la matrice */
	if (MatType(MatP3)!=Double)
	{
		Erreur("La matrice des coordonnees 3D doit etre Double");
		MatFree(&MatP3);
		return ERR_MAT_3D_TYPE;
	}
	/* Lecture de la matrice des points 2D */
	Matp2=MatReadAsc(argv[2]);
	if (Matp2==NULL)
	{
		Erreur("Probleme lors de la lecture de la matrice des coordonnees 2D");
		MatFree(&MatP3);
		return ERR_MAT_2D;
	}
	/* Verification du type de la matrice */
	if (MatType(Matp2)!=Double)
	{
		Erreur("La matrice des coordonnees 2D doit etre Double");
		MatFree(&MatP3);
		MatFree(&Matp2);
		return ERR_MAT_2D_TYPE;
	}

	/* Estimation de ma matrice de projection perspective */
	MatProjPers=Calibrage(MatP3,Matp2);
	if (MatProjPers==NULL)
	{
		Erreur("Probleme lors de l'estimation de la matrice de projection perspective");
		MatFree(&MatP3);
		MatFree(&Matp2);
		return ERR_MAT_OUT;
	}
	/* Ecriture de la matrice resultat */
	MatWriteAsc(MatProjPers,argv[3]);
	/* Liberation des matrices */
	MatFree(&MatP3);
    MatFree(&Matp2);
    MatFree(&MatProjPers);
    
	return OK;
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
    afficherMatrix(A);
    
    
    
    /* Création du vecteur colonne Val (valeur propres) 
    * et de la matrice des vecteurs propres Vect 
    */
    ata = Produit(Transposee(A), A);
    pVal = MatCAlloc(Double, MatNbCol(A), 1); 
    pVect = MatCAlloc(Double, MatNbCol(A),  MatNbCol(A));    
    
    SymEig(ata, &pVal, &pVect);
    
    afficherMatrix(pVect);
    
    
    
    
    /* Détermination de la valeur propre minimale */
    min_valP = min(pVal);
    printf("min=%d\n\n", min_valP);
    
    
    
    
    /* Création du vecteur propre M */
    m = MatCAlloc(Double, 3, 4); 
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
    /* Matrice avant normalisation */
    MatWriteAsc(m, "");

   
   
   
   /* Normalisation */
   norme = calcul_norme(M[2], 3);   
   for(i=0; i<MatNbRow(m); i++)
   {
    for(j=0; j<MatNbCol(m); j++)
    {
     M[i][j] = signe(M[MatNbRow(m)-1][MatNbCol(m)-1]) * M[i][j] / norme;
    }
   }
   /* Matrice après normalisation */
   MatWriteAsc(m, "");
   
   
   
   /* On libère la mémoire */
   MatFree(&A);
   MatFree(&ata);
   MatFree(&pVal);
   MatFree(&pVect);
   
   return m;
}

