function [D, L] = ALDLT(A)
    [n,n] = size(A);
    R = chol(A);

    D(i,i) = diag()
endfunction

exec('ARRT.sci', -1);

// SOLVE 

n = 5;
X = rand(n, n);
A = X*X';
// R = ARRT(A)
// norm(R'*R-X)
R = chol(A);
// norm(R'*R-A)

// FICHIER FONCTIONNEL
A = read('A.dat', 5,5)


// "ALGO" CM

// // D = diag(rii 2 ) et L = R T diag(rii )−1 .
// diagrii = zeros(n);
// diagri = diag(R.^R)
// for i = 1:n
//     diagrii(i, i) = diagri(i,1)
// end
// D = diagrii

// diagrii = zeros(n);
// diagri = diag(R)
// for i = 1:n
//     diagrii(i, i) = diagri(i,1)
// end
// L = R'*diagrii^-1



// ALGO WIKIPEDIA

function [D,L] = mychol(A)
    [n, n] = size(A);

    D = zeros(size(A));
    L = zeros(size(A));
    
    for j=1:n
        D(j,j) = A(j,j);
        v = 0;
        w = 0;
        for k=1:j-1
            v = v + A(j,k)^2 + D(k,k);
        end
        D(j,j) = D(j,j) - v;
    
        for i=j+1:n
            for k=i:j
                printf("%d, %d, %d\n", i, j, k);
                w = w + A(i,k)*A(j,k) * D(k,k);
            end
            L(i,j) = (A(i,j) - w) / D(j,j);
        end
    end
endfunction

// [D,L] = mychol(A)


function [D, L] = mychol2(A)
    [n, n] = size(A);

    L = zeros(n, n);
    d = zeros(n);
    v = zeros(n);

    for i = 1:n
        somme = 0;
        for j = 1 : i-1
            v(j) = L(i,j) * d(j);
            somme = somme + L(i,j) * v(j);
        end
        d(i) = A(i,i) - somme;
        for j = i+1:n
            somme = 0
            for k = 1:i-1
                somme = somme + L(j,k) * v(k)
            end
            L(j,i) = (A(j,i) - somme) / d(i);
        end
    end
    D = zeros(n,n);
    for j=1:n
        D(j,j) = d(j);
    end
    L = L + tril(eye(n,n));
endfunction

// [D,L] = mychol2(A)


// mychol2 REBORN
function [D, L] = mychol3(A)

    [n, n] = size(A);

    L = zeros(n, n);
    D = zeros(n,n);
    v = zeros(n);

    for i = 1:n
        somme = 0;
        for k = 1 : i-1
            v(k) = L(i,k) * D(k,k);
            // somme = somme + L(i,k) * v(k);
        end
        // D(i,i) = A(i,i) - somme;
        D(i,i) = A(i,i) - (L(i,1 : i-1) * v(1 : i-1));
        // D(i,i) = A(i,i) - L(1:i-1,i) * (L(i,1 : i-1) * D(1 : i-1,1 : i-1));
        for j = i+1:n
            somme = 0
            for k = 1:i-1
                somme = somme + L(j,k) * v(k)
            end
            L(j,i) = (A(j,i) - somme) / D(i,i);
        end
    end
    L = L + tril(eye(n,n));
endfunction

[D,L] = mychol2(A)
norm(L*D*L' - A)



