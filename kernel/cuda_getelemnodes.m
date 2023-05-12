function un = cuda_getelemnodes(un, u, np, ne, ncu, e1, e2)
% GPU implementation of un = u(:,e1:e2,:);

% [np, ne, ncu, e1, e2]
% size(u)
% size(un)

% nn = (e2-e1+1)*np;
% for i = 1:nn
%     for j=1:ncu
%         un(i+(j-1)*nn) = u(i+np*(e1-1)+(j-1)*np*ne);
%     end
% end
% 
% tm1=un(1:nn*ncu);
% tm2=reshape(u,[np ne ncu]);
% tm2=tm2(:,e1:e2,:);
% max(abs(tm1(:)-tm2(:)))

nn = (e2-e1+1)*np;
for i = 1:nn
    k = rem(i-1,np)+1;
    e = (i-k)/np+e1;  
    for j=1:ncu
        un(i+(j-1)*nn) = u(k+(j-1)*np+(e-1)*np*ncu);
    end
end

return;

np = 10;
ncu = 5;
ne = 100;
e1 = 21;
e2 = 50;
u = rand(np, ncu, ne);
un = zeros((e2-e1+1)*np,ncu);
un = cuda_getelemnodes(un, u, np, ne, ncu, e1, e2);
nn = (e2-e1+1)*np;
tm1=un(1:nn*ncu);
tm2=permute(u(:,:,e1:e2),[1 3 2]);
max(abs(tm1(:)-tm2(:)))

