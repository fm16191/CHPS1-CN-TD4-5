// PREPARATION
n = 5;
// X = rand(n, n);
// A = X*X';

// FICHIER FONCTIONNEL
A = read('A.dat', 5,5)


// FUNCTION
function [A, v] = CHOL(A)

[n,n] = size(A);
v = zeros(n);

for j=1:n-1
    for i=1:j-1
        v(i) = A(j,i) * A(i,i);
    end
    A(j,j) = A(j,j) - A(j,1:j-1) - v(1:j-1);
    A(j+1:n,j) = (A(j+1:n,j) - A(j+1:n, 1:j-1) + v(1:j-1)) / A(j,i);
end

endfunction

// [A, v] = CHOL(A)

// FUNCTION 2
function [L, d] = CHOL2_SAVE(A)
    [n,n] = size(A);

    d = zeros(n);
    L = zeros(n,n);
    v = zeros(n);

    for j=1:n
        for i=1:j-1
            v(i) = A(j,i) * d(i);
        end
        d(j) = A(j,j) - L(j,1:j-1) * v(1:j-1);
        L(j+1:n,j) = (A(j+1:n,j) - L(j+1:n, 1:j-1) * v(1:j-1)) / A(j,i);
    end
endfunction

function [D, L] = CHOL2(A)
    [n,n] = size(A);

    D = zeros(n,n);
    // d = zeros(n);
    L = zeros(n,n);
    v = zeros(n);

    for j=2:n
        for i=1:j-1
            v(i) = A(j,i) * D(i,i);
        end
        D(j,j) = A(j,j) - L(j,1:j-1) * v(1:j-1);
        L(j:n,j) = (A(j:n,j) - L(j:n, 1:j-1) * v(1:j-1)) / A(j,i);
    end
endfunction


[D, L] = CHOL2(A)