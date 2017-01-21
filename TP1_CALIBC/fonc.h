/*
*               Fichier contenant toutes les fonctions utiles 
*                       au calibrage d'une cam�ra
* 
*  par TERAITETIA Sara
*/

#ifndef __fonc_h_
#define __fonc_h_

#include "salade.h"
#include "limace.h"

/* Fonction d'affichage de matrice */
void afficherMatrix(Matrix Mat);

/* D�termine et renvoi la transpos�e d'une matrice */
Matrix Transposee(Matrix Mat);

/* Renvoi le produit de 2 matrices */
Matrix Produit(Matrix Mat1, Matrix Mat2);

/* Calcul la norme d'un vecteur */
double calcul_norme(double* vect, int len);

/* Fabrique la matrice A */
Matrix fabA(Matrix Mat3D, Matrix Mat2D);

/* D�termine l'indice du minimum d'une matrice */
int min(Matrix M);

/* Renvoi -1 ou 1 selon le signe du double mis en param�tre */
int signe(double m34);

#endif /* !__fonc_h_ */
