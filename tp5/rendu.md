
Réexpliquer l'exo de ce matin dans le code (9/12) (le travail préliminaire) : 

| --- | 
0 -> n+1 
pas : h




// TODO
> Besoin d'utilsiation de l'option de chargement de librairie `-lcblas` car [...]

p_env.c:(.text.startup+0x125) : référence indéfinie vers « cblas_dcopy »

> nécessité d'installation de cblas pour la dépendance `cblas.h`
For the best performance an optimized platform-specific CBLAS library should be used for -lcblas. The library must conform to the CBLAS standard. The ATLAS package provides a portable high-performance BLAS library with a CBLAS interface. It is free software and should be installed for any work requiring fast vector and matrix operations. The following command line will link with the ATLAS library and its CBLAS interface:
https://www.gnu.org/software/gsl/doc/html/usage.html

> Attention, potentiellement besoin de modifier des appels à LAPACK en conséquence



Architecture logicielle du projet 


# EX3

1. [Stockage](http://netlib.org/lapack/lug/node121.html)
1. [Stockage en bande](http://netlib.org/lapack/lug/node124.html)
Le stockage doit être en une seule dimension, de par la question 2, LAPACK_COL_MAJOR

[blas](http://www.netlib.org/blas/)

[blas routine reference](http://www.netlib.org/blas/blasqr.pdf)

COL major fonctionne niquel.

ROW major fonctionne, mais ne fonctionne dans Lapack ne marche pas (?) -> code bizarre
COL Major est efficace. => format le plus performant.


**Arrays are passed as pointers**, not as a pointer to pointers. All the LAPACKE routines that take **one or more 2D arrays as a pointer receive a single extra parameter of type int**. This argument must be equal to **either LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR** which are **defined in lapacke.h**, specifying whether the arrays are stored in **row-major or column-major order**. If a routine has multiple array inputs, they must all use the same ordering.

*Note that using row-major ordering may require more memory and time than column-major ordering, because the routine must transpose the row-major order to the column-major order required by the underlying LAPACK routine.*

Each 2D array argument in a FORTRAN LAPACK routine has an additional argument that specifies its leading dimension. For row-major 2D arrays, elements within a row are assumed to be contiguous and elements from one row to the next are assumed to be a leading dimension apart. For column-major 2D arrays, elements within a column are assumed to be contiguous and elements from one column to the next are assumed to be a leading dimension apart.

[netlib.org array arguments](http://www.netlib.org/lapack/lapacke.html#_array_arguments)

3. leading dimension ld

Dimensions principale
[Leading dimensions Usage](https://www.ibm.com/docs/en/essl/6.3?topic=matrices-how-leading-dimension-is-used)

La dimension dominaute c'est l'autre dimension que celle selon laquelle la matrice est stockée. Si la matrice est stockée par colonnes, alors la leading dimension est le nombre de lignes, donc n.

4. [dgbsv](https://www.ibm.com/docs/en/essl/6.2?topic=blaes-sgbsv-dgbsv-cgbsv-zgbsv-general-band-matrix-factorization-multiple-right-hand-side-solve) Implémente gauss avec pivots

5. Modifier pour stockage en priorité ligne et non colonne