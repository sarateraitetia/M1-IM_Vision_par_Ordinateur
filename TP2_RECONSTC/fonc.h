#ifndef __fonc_h_
#define __fonc_h_

#include "limace.h"
#include "salade.h"


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


#endif /* !__fonc_h_ */
