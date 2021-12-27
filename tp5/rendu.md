# **TP 5 : Résolution de l'équation de la chaleur en 1D stationaire**

Au cours de ce TP, nous cherchons à résoudre l'équation de la chaleur 1D stationnaire, avec implémentation d'algorithmes de résolution de systèmes linéaires.

# TP5 : Exercice 3 - Utilisation BLAS/LAPACK

Note : *De par des erreurs de compilations et de librairies manquantes sur le code de base, j'ai du installer `cblas`, et en accord avec la documentation [gnu](https://www.gnu.org/software/gsl/doc/html/usage.html), ajouter le flag de compilation `-lcblas`.*

> Attention, potentiellement besoin de modifier des appels à LAPACK en conséquence

1\. Pour utiliser Blas et Lapack, les matrices doivent être dans le format suivant: [¹](#references)
- tableau à deux dimensions
- stockée en paquet pour les matrices symétriques [²](#references)
- stockée en bandes pour les matrices bandes [³](#references)
- en vecteurs pour les matrices tridiagonales ou bidiagonales

2\. La constante LAPACK_COL_MAJOR spécifie que le tableau passé en argument est stocké en colonnes. [⁴](#references)

3\. La dimension principale est la taille de la première dimension de la matrice. Si on est en LAPACK_COL_MAJOR, par exemple, alors c'est le nombre de colonnes de cette matrice.

4\. dgbsv calcule la solution d'un système linéaire `A*X = B` pour des matrices de doubles stockées en général bande. Elle implémente une décomposition LU. [⁵](#references)

5\. La vérification de nos codes se fait dans le fichier modifié, après make et exécution de tp2poisson1D_direct en changeant le paramètre `row` dans `tp2_poisson_direct`, on obtient une erreur résiduelle relative d'ordre $10e^{-15}$ en COL et ROW major.
La vérification de nos fonctions se fait dans le fichier `tp2_poisson1D_direct.c`, par le calcul de l'erreur résiduelle entre le `RHS` calculé et la solution attendue, `EX_SOL`.

En row major, le code était déjà implémenté, et indiquait une erreur résiduelle de $10^{-15}$. Après implémentation de la fonction `set_GB_operator_rowMajor_poisson1D` dans `lib_poisson1D.c`, et re-exécution du programme avec `row = 1`, on obtient également une erreur résiduelle de l'ordre de $10^{-15}$. 

L'implémentation `write_GB_operator_colMajor_poisson1D` a aussi permie de confirmer que la matrice générée correspondait à celle générée en row. 

Note : *Dans l'appel à `LAPACKE_dgbsv`, la matrice de taille `lab*la` contient `kv` lignes/colonnes (en fonction de si row ou col major) supplémentaires à la matrice AB à initialiser à 0. Ces emplacements superflus sont utilisés comme vecteur de calcul par la fonction `LAPACKE_dgbtrf`, appelée par `LAPACKE_dgbsv`. Cette étape est nécessaire, pour le bon fonctionnement de la `dgbsv`.*


# TP5 : Exercice 4  - DGBMV

> DGBMV implémente une multiplication matrice vecteur sur un stockage Général Bande.
$y = \alpha A x + \beta y$

A été implémenté : 

```c
// Création d'un vecteur de test y pour vérifier l'utilisation de DGBMV
double *y;
y=(double *) malloc(sizeof(double)*la);
for (int i = 0 ; i < la; i++) y[i] = i+1;

// Utilisation
cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, 0, AB, la, y, 1, 0, y, 1);
// * Paramètres : 
// Entrée        : CblasColMajor - stockage de A en col major
// Entrée        : CblasNoTrans - A n'a pas besoin d'être transposée.
// Entrée        : la, la - m*n, c'est à dire nombre de lignes et colonnes
// Entrée        : kl - nombre de sousdiagonales (1 dans le cas d'une matrice tridiagonale)
// Entrée        : ku - nombre de surdiagonales (1 dans le cas d'une matrice tridiagonale)
// Entrée        : 0 - alpha de l'équation.
// Entrée        : AB - la matrice de l'équation
// Entrée        : la - dimension principale de A.
// Entrée        : y - vecteur x de l'équation.
// Entrée        : 1 - incX. pas d'itération du vecteur x.
// Entrée        : 0 - beta de l'équation.
// Entrée/Sortie : y - vecteur y de l'équation.
// Entrée        : 1 - incY. pas d'itération du vecteur y.

// Sauvegarder le résultat de y
write_vec(y, &la, "dgbmv_y_col_a.dat");
```

Afin de vérifier l'utilisation de `cblas_dgbmv`, on fait varier avec les coefficients $\alpha$ et $\beta$. On teste alors les équations suivantes : 
- $\alpha=0$ et $\beta=0$. 
  - Le résultat $y$ attendu est alors $y = 0$
  - Le résultat obtenu est en effet $0$.
- $\alpha=0$ et $\beta=1$. 
  - Le résultat $y$ attendu est alors $y = 0* A x + 1*y = y$
  - Le résultat obtenu est en effet $y$. En changeant la valeur de $\beta$ pour 2 par exemple, on obtient correctement $y = \beta y$.
- $\alpha=1$ et $\beta=0$. 
  - Le résultat $y$ attendu est alors $y = 1*A x + 0* y = Ax$
  - Malheuresement, le résultat obtenu est un vecteur nul.
  Sur le principe, la multiplication $AB*y$ est correcte sur toutes les valeurs sauf la première et dernière. $0$ est en effet attendu car $-1 * 1 + 2 * 1 + (-1) * 1 = 0$ mais en première et dernière positon sont attendues les valeurs $-1*1+2*1 = 1$, ce qui n'était malheuresement pas le cas. En changeant les valeurs initiales de y, un vecteur vide est tout de même renvoyé.
  Malgré le temps passé dessus, je n'ai pas réussi à comprendre pourquoi un tel résultat ou alors où était mon erreur.
- $\alpha=1$ et $\beta=1$.
  - Le résultat $y$ attendu est alors $y = 1*A x + 1* y = Ax+y$
  - Eest renvoyé un vecteur aléatoire, vecteur unitaire 1, dont les 5 premières valeurs ne correspondent à rien d'attendu. Cela n'est pas très étonnant, dans la mesure où le test précédent a échoué.

La même batterie de tests a aussi été effectuée en **row major**, et les mêmes résultats ont été obtenus, à savoir les deux premiers tests sont positifs, le troisième est un vecteur vide, et le 4ème un vecteur résutat incohérent.

> Bien que le choix de passer en argument une taille de matrice `M*N = la * la`, me semble pertinent étant donné qu'il s'agit d'une matrice stockée en général bande, je pense estimer mon erreur à cet endroit, même si après multiples tests, la configuration actuelle est la plus correcte à mes yeux d'un point de vue logique et testée.


# References
1. [Documentation Lapack @netlib.org - Schémas de stockage des matrices](http://netlib.org/lapack/lug/node121.html)
2. [Documentation Lapack @netlib.org - Packed storage](https://www.netlib.org/lapack/lug/node123.html)
3. [Documentation Lapack @netlib.org - Stockage en bande](http://netlib.org/lapack/lug/node124.html)
4. [constantes Lapack col/row major @netlib.org](http://www.netlib.org/lapack/lapacke.html#_array_arguments)
5. [Lapack's DGBSV @netlib.org](http://www.netlib.org/lapack/explore-html/d3/d49/group__double_g_bsolve_gafa35ce1d7865b80563bbed6317050ad7.html)




6. [Documentation Blas @netlib.org](http://www.netlib.org/blas/)
7. [Blas routine reference @netlib.org](http://www.netlib.org/blas/blasqr.pdf)