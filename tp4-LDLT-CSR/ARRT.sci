// CHOLESKY
// W=rand(5,5)+%i*rand(5,5);
// X=W*W';
// R=chol(X);
// norm(R'*R-X)

function [R, RT] = ARRT(A)
    [n,n] = size(A);
    R = zeros(n*n);

    for j = 1 : n
        for i = 1 : j - 1
            R(i, j) = A(i, j) - R(1 : (i - 1), i)' * R(1 : (i - 1), j);
            R(i, j) = R(i, j)/R(i, i);
            end
        R(j, j) = A(j, j) - R(:, j)' * R(:, j);
        R(j, j) = sqrt(R(j, j));
    end
endfunction

