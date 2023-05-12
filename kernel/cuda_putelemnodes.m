function u = cuda_putelemnodes(u, un, np, ne, ncu, e1, e2)
% GPU implementation of u(:,:,e1:e2) = un

% nn = (e2-e1+1)*np;
% for i = 1:nn
%     for j=1:ncu
%         u(i+np*(e1-1)+(j-1)*np*ne) = un(i+(j-1)*nn);
%     end
% end

nn = (e2-e1+1)*np;
for i = 1:nn
    k = rem(i-1,np)+1;
    e = (i-k)/np+e1;  
    for j=1:ncu
         u(k+(j-1)*np+(e-1)*np*ncu) = un(i+(j-1)*nn);
    end
end
