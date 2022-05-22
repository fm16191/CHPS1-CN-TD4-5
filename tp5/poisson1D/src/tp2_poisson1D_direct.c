/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"

int main(int argc,char *argv[]){
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv;
  int info;
  int NRHS;
  double T0, T1;
  double *RHS, *EX_SOL, *X;
  double *AB;
  double *y;

  double temp, relres;

// KV ? - ne marche pas si pas de KV à 1
// indice dans la doc de lapack
// lignes temporaires nécessaires pour dgbtrf
// lié à dgbtrf
// KL, KU dans la doc

// stockage general band
// Matrice A, tridiagonale, avec diagonale principale à 2, et les deux autres à -1.
// Matrice AB tridiagonale 2 / -1 : 
//  0 -1 ... -1
//  2  2 ...  2
// -1 -1 ...  0
// Stockage en bande : par row, on a donc [0, 2, -1].
// On ajoute un 0 de chaque côté [0, 0, 2, -1, 0] car nécessaire pour x raisons => kv l'indice est maintenant à 1

  NRHS=1;
  nbpoints=102;
  la=nbpoints-2; // la : taille de matrice
  T0=-5.0;
  T1=5.0;

  printf("--------- Poisson 1D ---------\n\n");

  // Initialisation
  RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);
  
  // FOR TESTS PURPOSES ONLY
  y=(double *) malloc(sizeof(double)*la);
  for (int i = 0 ; i < la; i++) y[i] = 1;
  // y=(double *) malloc(sizeof(double)*3);
  // for (int i = 0 ; i < 3; i++) y[i] = i+1;

  // Vecteur linéaire uniforme de ]0 à 1[
  set_grid_points_1D(X, &la);
  
  // Matrice vide avec premier à T0 et dernier à T1
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);

  // Vecteur linéaire uniforme de ]TC0 à TC1[
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  // Sauvegarde
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=1; // Deuxième indice
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

  // Matrice de sortie 
  AB = (double *) malloc(sizeof(double)*lab*la);

  // Résultat de LAPACKE_dgbsv
  info=0;

  /* working array for pivot used by LU Factorization */
  ipiv = (int *) calloc(la, sizeof(int));

  int row = 0; // choix de la méthode de Lapack

  if (row == 0){ // LAPACK_ROW_MAJOR
    printf("LAPACK_ROW_MAJOR\n");
    printf("lab : %d _ la : %d\n", lab, la);

    set_GB_operator_rowMajor_poisson1D(AB, lab, la);
    write_GB_operator_rowMajor_poisson1D(AB, lab, la, "AB_row.dat");

    info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR,la, kl, ku, NRHS, AB, la, ipiv, RHS, NRHS);


    set_GB_operator_rowMajor_poisson1D(AB, 3, la);
    write_vec(y, &la, "dgbmv_y_row_b.dat");
    // cblas_dgbmv(CblasRowMajor, CblasNoTrans, M, N, kl, ku, alpha, A, lda, X, incX, beta, Y, incY);
    cblas_dgbmv(CblasRowMajor, CblasNoTrans, la, la, kl, ku, 2, AB, la, y, 1, 1, y, 1);
    write_vec(y, &la, "dgbmv_y_row_a.dat");
  } 
  
  else { // LAPACK_COL_MAJOR
    printf("LAPACK_COL_MAJOR\n");
    set_GB_operator_colMajor_poisson1D(AB, 3, la, 0);
    write_GB_operator_colMajor_poisson1D(AB, lab, la, "AB_col.dat");
    // write_GB_operator_rowMajor_poisson1D(AB, &lab, &la, "AB_row.dat");

    // https://www.netlib.org/lapack/explore-html/de/ddd/lapacke_8h_a84fd8b92816876c9c6c20608a1362aa5.html
    // Ax = B
    // B : RHS
    // LDAB : leading dimension of AB
    // Ipiv : pivot de dimensions N
    // The pivot indices that define the permutation matrix P;
    // row i of the matrix was interchanged with row IPIV(i).

    // la : leading dimension of RHS
    // info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N, kl, ku, NRHS, AB, LDAB, ipiv, B, ldb);
    info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, la, kl, ku, NRHS, AB, lab, ipiv, RHS, la);


    printf("lab : %d _ la : %d\n", lab, la);

    // m = la = 100
    // n = lab = 3 <- LEADING DIMENSION (LDA)
    // Since we are in col major, 
    
    set_GB_operator_colMajor_poisson1D(AB, 3, la, 0);
    write_GB_operator_colMajor_poisson1D(AB, 3, la, "dgbmv_AB_col.dat");
    // cblas_dgbmv(CBLAS_LAYOUT layout,
    // CBLAS_TRANSPOSE TransA, const CBLAS_INDEX M, const CBLAS_INDEX N,
    // const CBLAS_INDEX KL, const CBLAS_INDEX KU, const double alpha,
    // const double* A, const CBLAS_INDEX lda, const double* X,
    // const CBLAS_INDEX incX, const double beta, double* Y, const CBLAS_INDEX incY);
    // cblas_dgbmv(CblasColMajor, CblasNoTrans, M, N, kl, ku, alpha, A, lda, X, incX, beta, Y, incY);

    printf("kl : %d _ ku : %d\n", kl, ku);

    write_vec(y, &la, "dgbmv_y_col_b.dat");
    cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, 1, AB, la, y, 1, 1, y, 1);
    write_vec(y, &la, "dgbmv_y_col_a.dat");
    // cblas_dgbmv(CblasRowMajor,
    // Computes alpha*A*x + beta*y or alpha*A'*x + beta*y depending on the value of TransA.
    
    // Fonctionne en colMajor, ||b - A*x||
  }    

  
  printf("\n INFO DGBSV = %d\n",info);

  // Ecriture solution RHS
  write_xy(RHS, X, &la, "SOL.dat");

  /* Relative residual */
  temp = cblas_ddot(la, RHS, 1, RHS,1);
  temp = sqrt(temp);
  cblas_daxpy(la, -1.0, RHS, 1, EX_SOL, 1);
  relres = cblas_ddot(la, EX_SOL, 1, EX_SOL,1);
  relres = sqrt(relres);
  relres = relres / temp;
  
  printf("\nThe relative residual error is relres = %e\n",relres);

  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  free(ipiv);

  printf("\n\n--------- End -----------\n");
  return 0;
}