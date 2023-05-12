function test_mQR()

positions = perms([1 2 3 4]);

size(positions)

%for posi = 1:size(positions,1)
for posi = 1:1
pos = positions(posi,:);

iw = pos(1);
ik = pos(2);
it = pos(3);
ie = pos(4);

label = {};

label{iw} = 'iw';
label{ik} = 'ik';
label{it} = 'it';
label{ie} = 'ie';

label

n(iw) = 128;
n(ik) = 128;
n(it) = 128;
n(ie) = 1024;

A = rand(n(iw),n(it),n(ie));
Q = A;
R = A;
Ainv = rand(n(iw),n(it),n(ie));

disp('MATLAB QR inside a loop')
tic
for i = 1:n(ie)
    [Q(:,:,i),R(:,:,i)] = qr(A(:,:,i));   
end
toc

stop

% disp('Our scalar QR inside a loop')
% tic
% for i = 1:n(ie)
%     [Q,R] = QRMGS(A(:,:,1));
% end
% toc

% Q1
% R1
% 
% Q
% R
% 
% max(Q*R - A(:,:,1))


mA = mArray.fromArray(n,A,[iw it ie]);

% size(mA)

tic
[mQ,mR] = mQR(n,mA,iw,it,ik,ie);
toc;

% Q(:,:,1)
% R(:,:,1)
% mArray.toArray(mQ(ie,1),[iw ik])
% mArray.toArray(mR(ie,1),[ik it])

end

function [Q,R] = mQR2(n,A,i,j,k,l)
Q = mArray.zeros(n,[i k l]); % TODO: Improve this initialization
R = mArray.zeros(n,[k j l]);

for lv = 1:n(l)
    [q r] = qr(mArray.toArray(A(l,lv),[i j]));
    
    Q(l,lv) = mArray.fromArray(n,q,[i k]);
    R(l,lv) = mArray.fromArray(n,r,[k j]);
end

function [Q,R] = mQR(n,A,i,j,k,l)
ni = n(i);

V = A;
Q = mArray.zeros(n,[i k l]); % TODO: Improve this initialization
R = mArray.zeros(n,[k j l]);

for iv = 1:ni   
   Viv = V(j,iv);
   Riviv = mNorm(Viv,i);   
   Qiv = Viv ./ Riviv;   
   R(k,iv,j,iv) = Riviv;
   
   jv = iv+1:ni;
   
   if numel(jv) > 0
       Vjv = V(j,jv);
       [Rivjv,ops,mops]= contract(Qiv,Vjv,i);
       Vjv = Vjv - Rivjv .* Qiv;
       R(k,iv,j,jv) = Rivjv;
       V(j,jv) = Vjv;
   end
   
   Q(k,iv) = Qiv;
end

% for iv = 1:ni   
%    R(k,iv,j,iv) = mNorm(V(j,iv),i);
%    Q(k,iv) = V(j,iv) ./ R(k,iv,j,iv);
%    
%    for jv = iv+1:ni
%        R(k,iv,j,jv) = contract(Q(k,iv),V(j,jv),i);
%        V(j,jv) = V(j,jv) - R(k,iv,j,jv).*Q(k,iv); 
%    end
% end

% for i = 1:n
%     R(i,i) = norm(V(:,i));
%     Q(:,i) = V(:,i)./R(i,i);
%     
%     for j = i+1:n
%         R(i,j) = Q(:,i)'*V(:,j);
%         V(:,j) = V(:,j) - R(i,j).*Q(:,i);
%     end
% end

function [norm] = mNorm(A,i)

norm = sqrt(contract(A,A,i));

function [Q,R] = QRMGS(A)

[n,n] = size(A);
V = A;
Q = zeros(n,n);
R = zeros(n,n);

for i = 1:n
    R(i,i) = norm(V(:,i));
    Q(:,i) = V(:,i)./R(i,i);
    
    for j = i+1:n
        R(i,j) = Q(:,i)'*V(:,j);
        V(:,j) = V(:,j) - R(i,j)*Q(:,i);
    end
end

% % function [] = test_mArray()
% 
% % Positions for the indices number 1,...,8
% % Changing this vector we change the positions of the indices 
% % in the whole of the code, without changing the code!
% pos = [1 2 3 4 5 6 7 8];
% 
% % Number of indices
% nIndices = numel(pos);
% 
% % Labeling the indices by positions (this is the additional level of indirection)
% % that allows to do the order-oblivious indexing in the rest of the code
% ie = pos(1); % element
% is = pos(2); % side
% ig = pos(3); % integration point
% it = pos(4); % trial function
% iw = pos(5); % weight / test function
% id = pos(6); % reference dimension
% ix = pos(7); % spatial dimension
% ic = pos(8); % solution component
% 
% % Generatig the vector n that stores the sizes of each index
% n = ones(1,nIndices);
% 
% % Putting the the sizes
% n(ie) = 5; % Five elements
% n(is) = 3; % 3 sides per element
% n(ig) = 28; % 28 integration points
% n(it) = 21; % 21 trial functions
% n(iw) = 21; % 21 weight functions
% n(id) = 2; % 2 variables for partial differentation in the reference element
% n(ix) = 2; % 2 spatial dimensions
% n(ic) = 4; % 4 components of the solution
% 
% % M = mArray(n,[ig it iw ic]);
% % J = mArray(n,[id ix ie ig]);
% 
% % Initialization of a Jacobian with ones depending on ix, id, ig, ie
% J = mArray.ones(n,[ix id ig ie]);
% 
% % The following lines are equivalent
% J
% J(ix,:,id,:)
% J(id,:,ix,:)
% J(ix,1:end,id,1:end); 
% J(id,1:end,ix,1:end); 
% 
% % Filling the Jacobian for all ig and ie (it works also for one element, or
% % one gauss point) !!!! This important because the computational code
% % of different types of elements could be the same !!!!
% J(ix,1,id,1) = 2;
% J(ix,2,id,2) = 2;
% J(ix,1,id,2) = 1;
% J(ix,2,id,1) = 1;
% 
% % Computing the determinant. Look inside mDet2x2
% detJ = mDet2x2(J,ix,id);
% 
% % Computing the inverse of the Jacobian. Look inside mInv2x2
% Jinv = mInv2x2(J,ix,id);
% 
% % Contracting the reference dimensions to test that Jinv is the inverse
% I = contract(J,Jinv,id);
% 
% % Should help to understand how to do fluxes
% F = mArray.ones(n,[ie ig ic]);
% G = F(ic,1) .* sin(F(ic,2)) .^ F(ic,3) - F(ic,4);
% 
% % Should help to understand how to do contractions
% C = mArray.ones(n,[iw it ic ie]);
% U = mArray.ones(n,[iw ie ic]);
% G = contract(C,U,[iw]);
% 
% % Converting an mArray to a MATLAB array with indices ie, ic, it and iw in
% % the prescribed order
% CC = mArray.toArray(C,[ie ic it iw]);
% 
% % Converting a MATLAB array to an mArray depending on indices ie, ic, it
% % and iw.
% CCC = mArray.fromArray(n,C,[ie ic it iw]);
% 
% n(ie)
% n(ic)
% n(it)
% n(iw)
% 
% size(C)
% size(U)
% size(G)
% 
% size(C)
% size(CC)
% size(CCC)
% 
