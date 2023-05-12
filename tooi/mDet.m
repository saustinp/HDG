function d = mDet(A,i,j)

sizeA = size(A);
m = sizeA(i);
n = sizeA(j);

if ( m == 2 && n == 2)
    d = mDet2x2(A,i,j);
elseif ( m == 3 && n == 3) 
    d = mDet3x3(A,i,j);
else
    disp('Not implemented for m and n different to 2');
end

function d = mDet2x2(A,i,j)

d = A(i,1,j,1) .* A(i,2,j,2) - A(i,1,j,2) .* A(i,2,j,1);

function d = mDet3x3(A,i,j)

d = A(i,1,j,1).*A(i,2,j,2).*A(i,3,j,3) - A(i,1,j,1).*A(i,3,j,2).*A(i,2,j,3)+ ...
    A(i,2,j,1).*A(i,3,j,2).*A(i,1,j,3) - A(i,2,j,1).*A(i,1,j,2).*A(i,3,j,3)+ ...
    A(i,3,j,1).*A(i,1,j,2).*A(i,2,j,3) - A(i,3,j,1).*A(i,2,j,2).*A(i,1,j,3);            


