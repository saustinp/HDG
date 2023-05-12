function [p,t] = cylinder_mesh

m = 50;
n = 30;


[Y,X] = meshgrid(0:1/n:1,0:1/m:1);

p = [X(:),Y(:)];


dlay = 1.0;
rat = 1.2;
d0 = 1e-4;
p(:,2) = p(:,2)*(n+1);
p(:,2) = d0*(rat.^p(:,2)-1)/(rat-1);
p(:,2) = p(:,2)*dlay/p(end,2);

t = [1, 2, m+3, m+2];
t = kron(t,ones(m,1)) + kron(ones(size(t)),(0:m-1)');
t = kron(t,ones(n,1)) + kron(ones(size(t)),(0:n-1)'*(m+1));
level = kron((n:-1:1)',ones(1,m)); level = level(:);

t2t = mkt2t(t,1);
figure(1);
meshplot(p,t);
grid on;

nt = size(t,1);
for ie = 1:nt
    if (level(ie) == 18),
       [p,t,t2t,level] = refine(ie,p,t,t2t,level);
    end
end 


nt = size(t,1);
for ie = 1:nt
    if (level(ie) == 22),
       [p,t,t2t,level] = refine(ie,p,t,t2t,level);
    end
end 


nt = size(t,1);
for ie = 1:nt
    if (level(ie) == 26),
       [p,t,t2t,level] = refine(ie,p,t,t2t,level);
    end
end

figure(2);
meshplot(p,t);
grid on;

pmap(:,1) = (1 + p(:,2)*4).*cos(p(:,1)*pi-0.5*pi);
pmap(:,2) = (1 + p(:,2)*4).*sin(p(:,1)*pi-0.5*pi);

figure(3);
meshplot(pmap,t);
grid on;