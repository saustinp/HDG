%function [] = Example()

%% 
% The following prefixes are used in this code:
% z:  order oblivious index (used for geometry)
% zi: related with test functions
% zj: related with trial functions
% zk: shared by test and trial (usually is going to be
% contracted)
% n:  number of entities of the type determined by the suffix (e.g. ne is
% number of elements)
%
% Meaning of suffixes:
% e: element
% f: face
% c: component, zic: equation component, zjc: unkwown component
% r: reference space dimension (master element space)
% x: physical space dimension
% s: shape function
% g: integration points

% Zeros are disabled indices
 zis =  1; zir =  2; zix =  0; zig = 0; zie =  0; zif =  0; zic =  0; % Related with test functions
 zjs =  9; zjr =  3; zjx =  0; zjg = 0; zje =  0; zjf =  0; zjc =  0; % Related with trial functions
 zks =  0; zkr =  0; zkx =  4; zkg = 0; zke =  0; zkf =  0; zkc =  0; % Shared by test and trial functions
  zs =  5;  zr =  6;  zx =  7;  zg = 8;  ze = 10;  zf =  0;  zc = 11; % Related with geometry mappings
  
% Labels for used indices
labels = {};
labels{zis} = 'zis'; labels{zir} = 'zir';
labels{zjs} = 'zjs'; labels{zjr} = 'zjr';
labels{zkx} = 'zkx';
labels{zs} = 'zs'; labels{zr} = 'zr'; labels{zx} = 'zx'; labels{zg} = 'zg'; labels{ze} = 'ze';

%%
% Creating some mock TOOIs to emulate computations later
% Examples of initizalizations
% Creating sizes:

ne = 100; % Number of elements
ns = 15;  % Number of shape functions
ng = 30;  % Number of integration points
nr = 2;   % Number of reference space dimensions (chi-space)
nx = 2;   % Number of physical space dimensions (x-space)

% Node coordinates initialization (mock)
% Order is irrelevant
% e.g. next line is equivalent to: X = mArray.ones([ns nx ne],[zs zx ze])

X = mArray.rand([ne ns nx],[ze zs zx]); 

% Integration weights (mock)
w = mArray.ones(ng,zg);

% Shape basis functions evaluated at integration points (mock)
Phi_chi = mArray.ones([ns ng],[zs zg]);

% Copying Phi_chi to create test and trial basis functions at chi
% Example of map function (maps indices to create a copy with different indices)
Phi_i_chi = map(Phi_chi,zs,zis); % Phi_i_chi depends on zis instead of zs
Phi_j_chi = map(Phi_chi,zs,zjs); % Phi_j_chi depends on zjs instead of zs

% Show indices for Phi_i_chi and Phi_j_chi
disp('Indices for Phi_i_chi');
indices(Phi_i_chi,labels)
disp('Indices for Phi_j_chi');
indices(Phi_j_chi,labels)

% Shape basis functions derivatives respect reference variables (r index)
% Derivatives are evaluated at integration points (mock)
DPhi_chi = mArray.ones([ns ng nr],[zs zg zr]);
DPhi_i_chi = map(DPhi_chi,[zs zr],[zis zir]); % Maps [zs zr] to [zis zir]
DPhi_j_chi = map(DPhi_chi,[zs zr],[zjs zjr]); % Maps [zs zr] to [zjs zjr]

% Show indices
% Example of indices function
disp('Indices for DPhi_i_chi');
indices(DPhi_i_chi,labels)
disp('Indices for DPhi_j_chi');
indices(DPhi_j_chi,labels)

%%
% Isoparametric mapping at integration points
% Example of contraction
X_chi = contract(X,Phi_chi,zs);
disp('Indices for X_chi');
indices(X_chi,labels)

%%
% Jacobian of the isoparametric mapping
J_chi = contract(X,DPhi_chi,zs);
disp('Indices for J_chi');
indices(J_chi,labels)

%%
% Determinant of the Jacobian of the isoparametric mapping
% Example of mDet
detJ_chi = mDet(J_chi,zx,zr);
disp('Indices for detJ_chi');
indices(detJ_chi,labels)

%%
% Inverse of the Jacobian of the isoparametric mapping
% Example of inverse function mInv
JInv_chi = mInv(J_chi,zx,zr,zr,zx);
disp('Indices for JInv_chi');
indices(JInv_chi,labels)

%% 
% Mapping the inverse of the Jacobian for test and trial use
JInv_i_chi = map(JInv_chi,[zx zr],[zkx zir]);
JInv_j_chi = map(JInv_chi,[zx zr],[zkx zjr]);
disp('Indices for JInv_i_chi');
indices(JInv_i_chi,labels)
disp('Indices for JInv_j_chi');
indices(JInv_j_chi,labels)

%% 
% Mass matrix for all the elements at the same time
M = contract(w .* Phi_i_chi .* Phi_j_chi,detJ_chi,zg);
disp('Indices for M');
indices(M,labels)

%%
% (uh,div v)
C = contract(w .* DPhi_i_chi .* Phi_j_chi,detJ_chi,[zkx zir zg]);
disp('Indices for C');
indices(C,labels)

%%
% Stiffness matrix
% Example of contraction of more than one index
K = contract(w .* DPhi_i_chi .* DPhi_j_chi, JInv_i_chi .* JInv_j_chi .* detJ_chi,[zkx zir zjr zg]);
disp('Indices for K');
indices(K,labels)

%% 
% Converting a TOOI to a standard MATLAB multi-array
% Example of toArray
M_ = mArray.toArray(M,[zis zjs ze]);

%% 
% Converting a standard MATLAB multi-array to a TOOI
% Example of fromArray
M2 = mArray.fromArray(M_,[zis zjs ze]);

% The following accesses are equivalent (access to the same data in the
% same order)
% Example of index accessing
J_chi(zx,1,zr,2); % Accessing to dx1 / dchi2
J_chi; 
J_chi(zx,:,zr,:);
J_chi(zr,:,zx,:);

% Setting values for all the elements
% Example of writing in a TOOI ( it puts 2 in zjs,1 zis,1 for all elements)
M(zjs,1,zis,1) = 2;

% Should help to understand how to do fluxes
% Example of element wise operations
nc = 5 % (5 components for F)
F = mArray.ones([ne ng nc],[ze zg zc]); % Mock function F
G = F(zc,1) .* sin(F(zc,2)) .^ F(zc,3) - F(zc,4);
disp('Indices for G');
indices(G,labels)




