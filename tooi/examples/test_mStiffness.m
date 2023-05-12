function [] = test_mArray()

% Positions for the indices number 1,...,8
% Changing this vector we change the positions of the indices 
% in the whole of the code, without changing the code!
qe = 1;
qf = 2;
qg = 3;
qt = 4;
qw = 5;
qr = 6;
qs = 7;
qx = 8;
qy = 9;
qc = 10;
qu = 11;

pos = [qe  qu qg  qr qx qy qs qc qf qw qt];

% Number of indices
nIndices = numel(pos);

% Labeling the indices by positions (this is the additional level of indirection)
% that allows to do the order-oblivious indexing in the rest of the code
ie = find(pos == qe); % element
jf = find(pos == qf); % side
ig = find(pos == qg); % integration point
it = find(pos == qt); % trial function
iw = find(pos == qw); % weight / test function
ir = find(pos == qr); % reference dimension for trial function
is = find(pos == qs); % reference dimension for weight / test function
ix = find(pos == qx); % spatial dimension for trial function
iy = find(pos == qy); % spatial dimension for weight / test function
ic = find(pos == qc); % solution component
iu = find(pos == qu); % shape function

% Generatig the vector n that stores the sizes of each index
n = ones(1,nIndices);

% 

p = 5;
np = (p+1)*(p+2)/2;
nc = 4;
nsd = 2;
ne = 8192;
nf = 3;

% Putting the the sizes
n(ie) = ne; % Five elements
n(jf) = nf; % 3 sides per element
n(ig) = np; % integration points
n(it) = np; % trial functions
n(iw) = np; % weight functions
n(ir) = nsd; % 2 variables for partial differentation in the reference element
n(is) = nsd; % 2 variables for partial differentation in the reference element
n(ix) = nsd; % 2 spatial dimensions
n(iy) = nsd; % 2 spatial dimensions
n(ic) = nc; % 4 components of the solution
n(iu) = np; % trial functions

labels = {};
labels{ie} = 'e';
labels{jf} = 'f';
labels{ig} = 'g';
labels{it} = 't';
labels{iw} = 'w';
labels{ir} = 'r';
labels{is} = 's';
labels{ix} = 'x';
labels{iy} = 'y';
labels{ic} = 'c';
labels{iu} = 'u';


%W = mArray.ones(n,[ig iw]);
%T = mArray.ones(n,[ig it]);
DW = mArray.ones(n,[ig iw is]); % Derivative of weight functions
DT = mArray.ones(n,[ig it ir]); % Derivative of trial functions
w = mArray.ones(n,ig); % Integration weights

X = mArray.ones(n,[ie ix iu]);
T = mArray.ones(n,[iu ig ir]);

C = mArray.ones(n,[ie ig ix iy]); % Bi-linear form (change point-by-point)

t = tic;

indices(X,labels)
indices(T,labels)
Jt = contract(X,T,iu);

Jw = mArray.ones(n,[ie ig is iy]);

w_DW_DT = w .* DW .* DT; % Does not depend on geometry

detJ = mDet2x2(Jt,ir,ix);
invDetJ = 1 ./ detJ;

invJDetJt = mInvJDetJ2x2(Jt,ir,ix);
invJDetJw = mInvJDetJ2x2(Jw,is,iy);

indices(invJDetJt,labels)
indices(C,labels)
invJDetJt_C = contract(invJDetJt,C,ix);
invJDetJt_C_invJDetJw = contract(invJDetJt_C,invJDetJw,iy);

invJt_C_invJw_detJ = invJDetJt_C_invJDetJw .* invDetJ;

%invJt_C_invJw_detJ = invJt_C_invJw_detJ(ig,1,ir,1,is,1);

t1 = tic;
[K,ops] = contract(invJt_C_invJw_detJ,w_DW_DT,[ig ir is]);
t1 = toc(t1)

GFLOPS = ops / 10^9 / t1

t = toc(t)

