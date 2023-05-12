function Ainv = mInvDetJ2x2(A,i,j)
% It is a symmetric function (swaps i by j and conversely)

Ainv = A;

Ainv(i,1,j,1) =  A(i,2,j,2);
Ainv(i,2,j,2) =  A(i,1,j,1);
Ainv(i,1,j,2) = -A(i,2,j,1);
Ainv(i,2,j,1) = -A(i,1,j,2);

% % Ainv = Ainv  d;