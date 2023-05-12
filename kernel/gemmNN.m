function C = gemmNN(C, A, B, alpha, beta, tid, m, n, k)
% GPU implementation of C = alpha*A*B + beta*C
% A(m, n), B(n,k), C(m,k)

i = floor(tid/m)+1;
j = floor(tid/k)+1;
i1 = (j-1)*m+i;
C(i1) = beta*C(i1);
for l=1:n     
    C(i1) = C(i1) + alpha*A((l-1)*m+i)*B((k-1)*n+l);
end

