function [D, E, F] = Jacobi(A)
    D = diag(A);
    E = tril(A, -1)
    F = triu(A, 1)
    B = zeros(A)


endfunction


// Résolution du système Ax = b

function [x]=Jacobi(A,b,x0,k)
//
// On veut résoudre le systeme Ax=b
// de manière récursive
// x0 condition initiale de l’algorithme de Jacobi
// k nombre d’itérations
//
    n=size(A,"c")
    x=x0;
    for p=1:k
        for i=1:n
            x(i)=(b(i)-A(i,1:n)*x+A(i,i)*x(i))/A(i,i)
        end
    end
endfunction


A=[264 110 98;110 105 21;98 21 54]
b=[12;8;-5]
x0=[0;0;0]
[x1]=Jacobi(A,b,x0,1000)

norm(A*x1-b)








function [x,e]=jacobi(A,b,n,epsilon)

     N= -(triu(A,1)+tril(A,-1))
     M= diag(diag(A))

    x=zeros(n,1);
    e=norm(A*x-b)


    while(e>epsilon)
     x = inv(M)*(N*x + b)
     e=norm(A*x-b)
    end

endfunction