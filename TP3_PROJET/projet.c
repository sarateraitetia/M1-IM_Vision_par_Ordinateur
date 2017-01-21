#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "fonc.h"

/*----------------------------------------------------------------------
* ----------------------- Programme Principal --------------------------
*-----------------------------------------------------------------------
*/
int main(int argc, char *argv[])
{     
    int x, y;
    
    
    /* Declaration des variables de la fonction main */
    Matrix MatP2_init=NULL, MatP2_final=NULL, MatProjPers=NULL;
    Image image_finale;
    unsigned char** Img_f;
    double **h;
    
  
    Image image = ImRead("photo.pgm");
    // Image image = ImRead("ma_photo.pgm");
    unsigned char** I = ImGetI(image);

    
	/* Création de la matrice des points initiaux */
	// fabrication de la matrice (u,v) -> (4,2)
	MatP2_init = MatReadAsc("p2d_init.mx"); 
    // MatP2_init = MatReadAsc("p2d_init_ma_photo.mx");
	
	/* Création de la matrice des points finaux */
    // dÃ©claration de la matrice (X,Y) -> (4,2)
	MatP2_final = MatReadAsc("p2d_final.mx"); 
    // MatP2_init = MatReadAsc("p2d_final_ma_photo.mx");

	/* Estimation de ma matrice de projection perspective */
	// fabrication de la matrice A
	// SymEigth
	// on rÃ©cupÃ¨re le vecteur normal associÃ© Ã  la valeur normale la plus petite
	MatProjPers = Calibrage(MatP2_final, MatP2_init);
	// -> donne le vecteur h
	
	
	/* Ecriture de la matrice resultat */
    //MatWriteAsc(MatProjPers, "p2d_calibMINE.mx");
	
	
	/* Calcul des pixels(u,v) résultats */
	image_finale = ImAlloc(GrayLevel, ImNbRow(image), ImNbCol(image));
	Img_f = ImGetI(image_finale);
	
	h = MatGetDouble(MatProjPers);
	
	int u, v;
    for(x=0; x<ImNbRow(image_finale); x++)
    {
      for(y=0; y<ImNbCol(image_finale); y++)
	  {
             u = ceil((h[0][0] * x + h[0][1] * y + h[0][2]) / (h[2][0] * x + h[2][1] * y + h[2][2])); 
             v = ceil((h[1][0] * x + h[1][1] * y + h[1][2]) / (h[2][0] * x + h[2][1] * y + h[2][2])); 
             
             Img_f[x][y] = I[u][v];
      } 
    }
    ImWrite(image_finale, "test.ppm");
	// ImWrite(image_finale, "ma_photo_res.ppm");
	
	/* Liberation des matrices */
	MatFree(&MatP2_init);
    MatFree(&MatP2_final);
    MatFree(&MatProjPers);
    /* Liberation de l'image */
    ImFree(&image);
	
    system("pause");
	return 0;
}
