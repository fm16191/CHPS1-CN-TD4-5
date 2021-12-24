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

// function [D,L] = myLDLT(A)
//     [n, n] = size(A);

//     D = zeros(size(A));
//     // L = zeros(size(A));
//     L = eye(n, n);
    
//     for j=1:n
//         D(j,j) = A(j,j);
//         v = 0;
//         somme = 0;
//         for k=1:j-1
//             v = v + A(j,k)^2 + D(k,k);
//         end
//         D(j,j) = D(j,j) - v;
    
//         for i=j+1:n
//             for k=i:j
//                 // printf("%d, %d, %d\n", i, j, k);
//                 somme = somme + A(i,k)*A(j,k) * D(k,k);
//             end
//             L(i,j) = (A(i,j) - somme) / D(j,j);
//         end
//     end
// endfunction


function [L, D] = myLDLT3b(A)
    n = size(A, "r");

    L = eye(n, n);
    d = zeros(n);
    v = zeros(n);

    for i = 1:n
        u = 0;
        for j = 1 : i-1
            v(j) = L(i,j) * d(j);
            u = u + L(i,j) * v(j);
        end
        d(i) = A(i,i) - u;
        for j = i+1:n
            w = 0
            for k = 1:i-1
                w = w + L(j,k) * v(k)
            end
            L(j,i) = (A(j,i) - w) / d(i);
        end
    end
    D = zeros(n,n);
    for j=1:n
        D(j,j) = d(j);
    end
endfunction

function [L, D] = myLDLT1b(A)
    n = size(A, "r");
    for k = 1 : n-1
        A(k+1:n,k) = A(k+1:n, k) / A(k,k)
        A(k+1:n, k+1 : n) = A(k+1:n, k+1 : n) - A(k+1:n, k)*A(k,k+1 : n)
    end

    L = tril(A, -1);
    for i = 1:n
        D(i, i) = A(i, i);
        L(i, i) = 1;
    end
endfunction
