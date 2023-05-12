function [  ] = david(  )


% positions = perms([1 2 3 4]);
% 
% size(positions)
% 
% for posi = 1:size(positions,1)
% %for posi = 1:1
% pos = positions(posi,:);

pos = [4 3 2 1];

iw = pos(1);
ik = pos(2);
it = pos(3);
ie = pos(4);

labels{iw} = 'w';
labels{ik} = 'k';
labels{it} = 't';
labels{ie} = 'e';

% label = {};
% 
% label{iw} = 'iw';
% label{ik} = 'ik';
% label{it} = 'it';
% label{ie} = 'ie';
% 
% label

n(iw) = 16;
n(ik) = 16;
n(it) = 16;
n(ie) = 65536;

A = rand(n(iw),n(it),n(ie));

[Q,R] = QRMGS(A(:,:,1));


A = mArray.rand(n,[iw it ie]);

tic
[Q,R] = tooiQR(A,n,iw,it,ik,ie);
toc

B = contract(Q,R,ik);

size(A-B)
sum(A - B,[iw it ie])

indices(A,labels)
indices(Q,labels)
indices(R,labels)

end

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

end

function [value] = tooiNorm(A,k)
value = sqrt(sum(A.*A,k));
end

function [Q,R] = tooiQR(A,pn,pi,pj,pk,pl)

n = pn(pi);
V = A;

Q = mArray.zeros(pn,[pi pk pl]);
R = mArray.zeros(pn,[pk pj pl]);

for i = 1:n
    R(pk,i,pj,i) = tooiNorm(V(pi,:,pj,i),pi);
    Q(pi,:,pk,i)  = V(pi,:,pj,i) ./ R(pk,i,pj,i);
    
    for j = i+1:n
        R(pk,i,pj,j) = contract(Q(pi,:,pk,i),V(pi,:,pj,j),pi);
        V(pi,:,pj,j) = V(pi,:,pj,j) - R(pk,i,pj,j).*Q(pi,:,pk,i);
    end 
end

end