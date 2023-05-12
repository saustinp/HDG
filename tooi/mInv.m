function Ainv = mInv(A,ii1,jj1,jj2,ii2) % \warning ii1, jj1, jj2, ii2 components of multi-index should be ordered in a compatible way
%%
% Computes the inverse of the matrices determined by i1 and j.
% Note that $A$ maps from i1 to j, and the inverse $A^{-1}$ maps from j to
% i2.

sizeA = size(A);
m = prod(sizeA(ii1));
n = prod(sizeA(jj1));

if ( m == 2 && n == 2)
    Ainv = mInv2x2(A,ii1,jj1,jj2,ii2); % \warning Inverse 2x2 doe s not allow multi-indices
elseif ( m == 3 && n == 3)
    Ainv = mInv3x3(A,ii1,jj1,jj2,ii2); % \warning Inverse 3x3 does not allow multi-indices    
else
    Ainv = mInvNxN(A,ii1,jj1,jj2,ii2);
end

function Ainv = mInv2x2(A,i1,j1,j2,i2) 
% It is a symmetric function (swaps i by j and conversely)

d = mDet2x2(A,i1,j1);

ll = findL(A,i1,j1,j2,i2);
sizesA = size(A);
ni2 = sizesA(i2);
nj2 = sizesA(j2);
nll = sizesA(ll);

Ainv = mArray.zeros([ni2 nj2 nll],[i2 j2 ll]);

Ainv(i2,1,j2,1) =  A(i1,2,j1,2) ./ d;
Ainv(i2,2,j2,2) =  A(i1,1,j1,1) ./ d;
Ainv(i2,1,j2,2) = -A(i1,2,j1,1) ./ d;
Ainv(i2,2,j2,1) = -A(i1,1,j1,2) ./ d;

% % Ainv = Ainv  d;

function Ainv = mInv3x3(A,i1,j1,j2,i2) 
% It is a symmetric function (swaps i by j and conversely)

d = mDet3x3(A,i1,j1);

ll = findL(A,i1,j1,j2,i2);
sizesA = size(A);
ni2 = sizesA(i2);
nj2 = sizesA(j2);
nll = sizesA(ll);

Ainv = mArray.zeros([ni2 nj2 nll],[i2 j2 ll]);

Ainv(i2,1,j2,1) = -(A(i1,2,j1,3).*A(i1,3,j1,2) - A(i1,2,j1,2).*A(i1,3,j1,3))./d;
Ainv(i2,1,j2,2) = -(A(i1,2,j1,1).*A(i1,3,j1,3) - A(i1,2,j1,3).*A(i1,3,j1,1))./d;
Ainv(i2,1,j2,3) = -(A(i1,2,j1,2).*A(i1,3,j1,1) - A(i1,2,j1,1).*A(i1,3,j1,2))./d;
Ainv(i2,2,j2,1) = -(A(i1,1,j1,2).*A(i1,3,j1,3) - A(i1,1,j1,3).*A(i1,3,j1,2))./d;
Ainv(i2,2,j2,2) = -(A(i1,1,j1,3).*A(i1,3,j1,1) - A(i1,1,j1,1).*A(i1,3,j1,3))./d;
Ainv(i2,2,j2,3) = -(A(i1,1,j1,1).*A(i1,3,j1,2) - A(i1,1,j1,2).*A(i1,3,j1,1))./d;
Ainv(i2,3,j2,1) = -(A(i1,1,j1,3).*A(i1,2,j1,2) - A(i1,1,j1,2).*A(i1,2,j1,3))./d;
Ainv(i2,3,j2,2) = -(A(i1,1,j1,1).*A(i1,2,j1,3) - A(i1,1,j1,3).*A(i1,2,j1,1))./d;
Ainv(i2,3,j2,3) = -(A(i1,1,j1,2).*A(i1,2,j1,1) - A(i1,1,j1,1).*A(i1,2,j1,2))./d;

function ll = findL(A,i1,j1,j2,i2)
maxIndex = max([ndims(getData(A)),i1,j1,j2,i2]);

r = ones(1,maxIndex);
r(i1) = 0;
r(j1) = 0;
r(i2) = 0;
r(j2) = 0;
indices = 1:maxIndex;
ll = indices(r > 0);

function Ainv = mInvNxN(A,ii1,jj1,jj2,ii2)
ll = findL(A,ii1,jj1,jj2,ii2);

A1 = mArray.toArray(A,[ii1 jj1 ll]);

numIi1 = numel(ii1);
numJj1 = numel(jj2);
numLl  = numel(ll);

numIndices = numIi1 + numJj1 + numLl;

sizesA1 = ones(1,numIndices);
sizesA1(1:ndims(A1)) = size(A1);

nii1 = sizesA1( 1:numIi1 );
njj1 = sizesA1( (numIi1+1):(numIi1+numJj1) );
nll  = sizesA1( (numIi1+numJj1+1):(numIi1+numJj1+numLl) );

ni1 = prod(nii1);
nj1 = prod(njj1);
nl  = prod(nll);

A2 = reshape(A1,[ni1 nj1 nl]);
A2inv = A2; 

for l = 1:nl
    A2inv(:,:,l) = inv(A2(:,:,l));
end

A1inv = reshape(A2inv,[nii1 njj1 nll]);
Ainv = mArray.fromArray(A1inv,[jj2 ii2 ll]);

% function Ainv = mInvNxNold(A,i1,j1,j2,i2) 
% ll = findL(A,i1,j1,j2,i2);
% 
% A1 = mArray.toArray(A,[i1 j1 ll]);
% 
% sizesA1 = size(A1);
% 
% ni = sizesA1(1);
% nj = sizesA1(2);
% nll = sizesA1(3:ndims(A1));
% nl = prod(sizesA1(3:ndims(A1)));
% 
% A2 = reshape(A1,[ni nj nl]);
% A2inv = A2; 
% 
% for l = 1:nl
%     A2inv(:,:,l) = inv(A2(:,:,l));
% end
% 
% A1inv = reshape(A2inv,[ni nj nll]);
% Ainv = mArray.fromArray(A1inv,[j2 i2 ll]);
