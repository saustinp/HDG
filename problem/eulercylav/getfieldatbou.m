function [x,y,u] = getfieldatbou(mesh,u)

x = mesh.dgnodes(:,1,:);
y = mesh.dgnodes(:,2,:);
ind = find(abs(x(:).^2 + y(:).^2)<1+1e-10);
[y,jj] = sort(y(ind));
x = x(ind(jj));
u = u(ind(jj));
figure(1); clf; plot(y(:), u(:));


% % p = mesh.porder;
% % x = reshape(x,p+1,[]);
% % y = reshape(y,p+1,[]);
% % 
% % z = linspace(0,1,100)';
% % nfs = mkshape(p,mesh.plocfc,z,1);
% % fs = nfs(:,:,1)';
% % x = fs*x;
% % y = fs*y;
% % x = x(:);
% % y = y(:);




