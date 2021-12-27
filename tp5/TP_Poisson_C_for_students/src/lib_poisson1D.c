/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include <string.h>

void set_GB_operator_rowMajor_poisson1D(double* AB, int lab, int la) {
  // TODO
  // Faire les variables kv, lab, la plutôt que des pointeurs (?)
  // printf("lab : %d _ la : %d\n", lab, la);
  int kv = lab - 3;
  int ii, jj, kk;
  // Si l'index dans la matrice différent de 0, on a kv colonnes à initialiser à 0.
  for (jj = 0; jj < kv; jj++) {
    // Pour chaque ligne index < kv, on initialise à 0
    for (ii = 0;ii < la; ii++) {
      AB[kv * lab + ii] = 0.0;
    }
  }

  // Pour chaque ligne index < kv, on initialise à 0
  for (ii = 0; ii < la; ii++) {
    AB[(kv)*la + ii] = -1.0;
    AB[(kv + 1) * la + ii] = 2.0;
    AB[(kv + 2) * la + ii] = -1;
  }

  AB[0] = 0.0;
  if (kv == 1) { AB[kv * la] = 0; }
  // Remplacer par -> 
  // AB[kv * la] = 0;

  AB[(lab) * (la)-1] = 0.0;
}

// Initialisation de la matrice
void set_GB_operator_colMajor_poisson1D(double* AB, int lab, int la, int kv) {
  // Faire les variables kv, lab, la plutôt que des pointeurs (?)
  // printf("lab : %d _ la : %d\n", lab, la);
  int ii, jj, kk;
  for (jj = 0;jj < (la);jj++) {
    kk = jj * (lab);
    // kk indentation dans la matrice kk = (jj) nombre de colonnes * (lab) taille colonne

    // Si l'index dans la matrice différent de 0, on a kv colonnes à initialiser à 0.
    if (kv >= 0) {
      for (ii = 0;ii < kv;ii++) { // Pour chaque ligne index < kv, on initialise à 0
        AB[kk + ii] = 0.0;
      }
    }

    AB[kk + kv] = -1.0;
    AB[kk + kv + 1] = 2.0;
    AB[kk + kv + 2] = -1.0;
  }
  AB[0] = 0.0;
  if (kv == 1) { AB[1] = 0; }

  AB[(lab) * (la)-1] = 0.0;

  // for (jj = 0;jj < (la);jj++) {
  //   printf("=> %lf, %lf, %lf, %lf\n", AB[jj * lab], AB[jj * lab + 1], AB[jj * lab + 2], AB[jj * lab + 3]);
  // }
}


void set_GB_operator_colMajor_poisson1D_Id(double* AB, int* lab, int* la, int* kv) {
  int ii, jj, kk;
  for (jj = 0;jj < (*la);jj++) {
    kk = jj * (*lab);
    if (*kv >= 0) {
      for (ii = 0;ii < *kv;ii++) {
        AB[kk + ii] = 0.0;
      }
    }
    AB[kk + *kv] = 0.0;
    AB[kk + *kv + 1] = 1.0;
    AB[kk + *kv + 2] = 0.0;
  }
  AB[1] = 0.0;
  AB[(*lab) * (*la) - 1] = 0.0;
}

// Matrice vide de taille la avec premier à BC0 et dernier à BC1
void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1) {
  int jj;
  RHS[0] = *BC0;
  RHS[(*la) - 1] = *BC1;
  for (jj = 1;jj < (*la) - 1;jj++) {
    RHS[jj] = 0.0;
  }
}

// Vecteur linéaire uniforme EX_SOL de taille la de ]BC0 à BC1[
void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1) {
  int jj;
  double h, DELTA_T;
  // différence entre TC1 et TC0
  DELTA_T = (*BC1) - (*BC0);
  for (jj = 0;jj < (*la);jj++) {
    EX_SOL[jj] = (*BC0) + X[jj] * DELTA_T;
  }
}

// Vecteur linéaire uniforme de taille la ]0 à 1[
void set_grid_points_1D(double* x, int* la) {
  int jj;
  double h;
  h = 1.0 / (1.0 * ((*la) + 1));
  // pas de 1/taille+1
  for (jj = 0;jj < (*la);jj++) {
    x[jj] = (jj + 1) * h;
  }
}


/// Sauvegarde des données

// Sauvegarde du tableau AB, de taille lab*la en ROW major
void write_GB_operator_rowMajor_poisson1D(double* AB, int lab, int la, char* filename) {
  FILE* file;
  int ii, jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL) {
    for (ii = 0;ii < (lab);ii++) {
      for (jj = 0;jj < (la);jj++) {
        fprintf(file, "%lf\t", AB[ii * (la)+jj]);
      }
      fprintf(file, "\n");
    }
    fclose(file);
  }
  else {
    perror(filename);
  }
}

// Sauvegarde du tableau AB, de taille lab*la en COL major
void write_GB_operator_colMajor_poisson1D(double* AB, int lab, int la, char* filename) {
  //TODO
  FILE* file;
  int ii, jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL) {
    for (jj = 0;jj < (la);jj++) {
      for (ii = 0;ii < (lab);ii++) {
        fprintf(file, "%lf\t", AB[jj * (lab)+ii]);
      }
      fprintf(file, "\n");
    }
    fclose(file);
  }
  else {
    perror(filename);
  }
}

// Sauvegarde du vecteur vec, de taille la
void write_vec(double* vec, int* la, char* filename) {
  int jj;
  FILE* file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL) {
    for (jj = 0;jj < (*la);jj++) {
      fprintf(file, "%lf\n", vec[jj]);
    }
    fclose(file);
  }
  else {
    perror(filename);
  }
}

// Sauvegarde du vecteur vec, de taille la avec son identifiant x
void write_xy(double* vec, double* x, int* la, char* filename) {
  int jj;
  FILE* file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL) {
    for (jj = 0;jj < (*la);jj++) {
      fprintf(file, "%lf\t%lf\n", x[jj], vec[jj]);
    }
    fclose(file);
  }
  else {
    perror(filename);
  }
}



void eig_poisson1D(double* eigval, int* la) {
  int ii;
  double scal;
  for (ii = 0; ii < *la; ii++) {
    scal = (1.0 * ii + 1.0) * M_PI_2 * (1.0 / (*la + 1));
    eigval[ii] = sin(scal);
    eigval[ii] = 4 * eigval[ii] * eigval[ii];
  }
}

double eigmax_poisson1D(int* la) {
  double eigmax;
  eigmax = sin(*la * M_PI_2 * (1.0 / (*la + 1)));
  eigmax = 4 * eigmax * eigmax;
  return eigmax;
}

double eigmin_poisson1D(int* la) {
  double eigmin;
  eigmin = sin(M_PI_2 * (1.0 / (*la + 1)));
  eigmin = 4 * eigmin * eigmin;
  return eigmin;
}

double richardson_alpha_opt(int* la) {
  //TODO
}

void richardson_alpha(double* AB, double* RHS, double* X, double* alpha_rich, int* lab, int* la, int* ku, int* kl, double* tol, int* maxit) {
  //TODO
}
