function [x,y] = getfieldaty(mesh,u)

x = mesh.dgnodes(:,1,:);
y = mesh.dgnodes(:,2,:);
ind = find(abs(y)<2e-3);
[x,jj] = sort(x(ind));
y = u(ind(jj));
x = x(1:2:end);
y = 0.5*(y(1:2:end) + y(2:2:end));

p = mesh.porder;
x = reshape(x,p+1,[]);
y = reshape(y,p+1,[]);

z = linspace(0,1,100)';
nfs = mkshape(p,mesh.plocfc,z,1);
fs = nfs(:,:,1)';
x = fs*x;
y = fs*y;
x = x(:);
y = y(:);

figure(1); clf; plot(x(:), y(:));


