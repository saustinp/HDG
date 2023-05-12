% function [] = test_mArray()

% Positions for the indices number 1,...,8
% Changing this vector we change the positions of the indices 
% in the whole of the code, without changing the code!
pos = [1 2 3 4 5 6 7 8];

% Number of indices
nIndices = numel(pos);

% Labeling the indices by positions (this is the additional level of indirection)
% that allows to do the order-oblivious indexing in the rest of the code
ie = pos(1); % element
is = pos(2); % side
ig = pos(3); % integration point
it = pos(4); % trial function
iw = pos(5); % weight / test function
id = pos(6); % reference dimension
ix = pos(7); % spatial dimension
ic = pos(8); % solution component

% Generatig the vector n that stores the sizes of each index
n = ones(1,nIndices);

% Putting the the sizes
n(ie) = 5; % Five elements
n(is) = 3; % 3 sides per element
n(ig) = 28; % 28 integration points
n(it) = 21; % 21 trial functions
n(iw) = 21; % 21 weight functions
n(id) = 2; % 2 variables for partial differentation in the reference element
n(ix) = 2; % 2 spatial dimensions
n(ic) = 4; % 4 components of the solution

% M = mArray(n,[ig it iw ic]);
% J = mArray(n,[id ix ie ig]);

% Initialization of a Jacobian with ones depending on ix, id, ig, ie
J = mArray.ones(n,[ix id ig ie]);

% The following lines are equivalent
J
J(ix,:,id,:)
J(id,:,ix,:)
J(ix,1:end,id,1:end); 
J(id,1:end,ix,1:end); 

% Filling the Jacobian for all ig and ie (it works also for one element, or
% one gauss point) !!!! This important because the computational code
% of different types of elements could be the same !!!!
J(ix,1,id,1) = 2;
J(ix,2,id,2) = 2;
J(ix,1,id,2) = 1;
J(ix,2,id,1) = 1;

% Computing the determinant. Look inside mDet2x2
detJ = mDet2x2(J,ix,id);

% Computing the inverse of the Jacobian. Look inside mInv2x2
Jinv = mInv2x2(J,ix,id);

% Contracting the reference dimensions to test that Jinv is the inverse
I = contract(J,Jinv,id);

% Should help to understand how to do fluxes
F = mArray.ones(n,[ie ig ic]);
G = F(ic,1) .* sin(F(ic,2)) .^ F(ic,3) - F(ic,4);

% Should help to understand how to do contractions
C = mArray.ones(n,[iw it ic ie]);
U = mArray.ones(n,[iw ie ic]);
G = contract(C,U,[iw]);

% Converting an mArray to a MATLAB array with indices ie, ic, it and iw in
% the prescribed order
CC = mArray.toArray(C,[ie ic it iw]);

% Converting a MATLAB array to an mArray depending on indices ie, ic, it
% and iw.
CCC = mArray.fromArray(n,C,[ie ic it iw]);

n(ie)
n(ic)
n(it)
n(iw)

size(C)
size(U)
size(G)

size(C)
size(CC)
size(CCC)

