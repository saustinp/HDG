function Ainv = mInv2x2(A,i,j) % \todo{Change signature A,i1,j,i2}
% It is a symmetric function (swaps i by j and conversely)

d = mDet2x2(A,i,j);

Ainv = A;

Ainv(i,1,j,1) =  A(i,2,j,2) ./ d;
Ainv(i,2,j,2) =  A(i,1,j,1) ./ d;
Ainv(i,1,j,2) = -A(i,2,j,1) ./ d;
Ainv(i,2,j,1) = -A(i,1,j,2) ./ d;

% % Ainv = Ainv  d;