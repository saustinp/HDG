function C = gemmTN(C, A, B, alpha, beta, tid, m, n, k)
% GPU implementation of C = alpha*A^T*B + beta*C
% A(n, m), B(n,k), C(m,k)

i = floor(tid/m)+1;
j = floor(tid/k)+1;
i1 = (j-1)*m+i;
C(i1) = beta*C(i1);
for l=1:n     
    C(i1) = C(i1) + alpha*A((i-1)*n+l)*B((k-1)*n+l);
end

