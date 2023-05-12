function d = mDet2x2(A,i,j)

d = A(i,1,j,1) .* A(i,2,j,2) - A(i,1,j,2) .* A(i,2,j,1);
