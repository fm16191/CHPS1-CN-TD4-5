// calculer le produit d'une matrice creuse A de m lignes et n colonnes
// avec un vecteur de longueur n

function [sp] = prodCSR(A, v)
    [m, n] = size(A);
    n2 = size(v, "r");
    if n ~= n2 then error("Matrix and vector must be same size")
; end
    sp = A * v
endfunction

// Real value : sp = A*v

function [AX, AI, AJ] = csmtCSR(A)
    [m,n] = size(A);
    pos = 1;
    AI = zeros(1,m+1);
    AJ = zeros(1,0);
    AX = zeros(1,0)
    for i=1:m
        for j=1:n
            if A(i,j) ~= 0 then
                AX = [AX, A(i,j)];
                AJ = [AJ, j];
                pos = pos + 1
            end
        end
        AI(i+1) = pos-1;
    end
    AX = full(AX)
endfunction


function [Av] = mCSRv(AX, AI, AJ, v)
    n = size(AI, 'c');
    Av = zeros(1,n-1);
    for i = 1:n-1
        jstart = AI(i);
        jend = AI(i+1);
        for j =jstart+1:jend
            Av(i) = Av(i) + AX(j) * v(AJ(j))
        end
    end
endfunction




// A = [0, 2, 0, 4, 0; 6, 0, 8, 0, 10; 11, 0, 0, 14, 0];
// [AX, AI, AJ] = csmtCSR(A)

// sp = sprand(3,5,0.2)
// [AX, AI, AJ] = csmtCSR(sp)


// v = [1;1;1;1;1]
// v = grand(5,1, "bin", 1,0.5)

// [av] = mCSRv(AX, AI, AJ, v)

// av == (A*v)'
