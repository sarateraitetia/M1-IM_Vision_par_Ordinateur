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


void afficherMatrix(Matrix Mat);

/* Produit renvoi le produit de => OK !!
 * 2 matrices
 */
Matrix Produit(Matrix Mat1, Matrix Mat2);

/* Transposee renvoi la matrice transposée => OK !!
 * de la matrice mis en paramètre
 */
Matrix Transposee(Matrix Mat);

/* renvoi l'indice du minimum 
* d'une matrice M
*/
int min(Matrix M);

int signe(double m34);

/* Norme vectorielle 
 *
 */
double calcul_norme(double* vect, int len);

/* Fabrication de la matrice A
 * 
 */
Matrix fabA(Matrix Mat3D, Matrix Mat2D);

/* Calibrage 3D : estimation au sens des moindres carres totaux de la matrice
 * de projection perspective à partir de correspondances 3D <--> 2D
 */
Matrix Calibrage(Matrix MatP3, Matrix Matp2);


#endif
