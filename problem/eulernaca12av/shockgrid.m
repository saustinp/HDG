function [p,t] = shockgrid(x1, y1, ds, ns, parity)

m = length(ds);
n = length(x1);
sx = zeros(n, m);
sy = 0*sx;
for i = 1:m
    sx(:,i) = x1(:) + ds(i)*ns(1);
    sy(:,i) = y1(:) + ds(i)*ns(2);
end
sx=sx';
sy=sy';
p = [sx(:) sy(:)];

% t = [1 2 m+2 m+1];
% t = kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m);
% t = kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');
% fix p
for i = 1:size(p,1)
    p0 = p(i,:);
    d1 = p0(2).^2-(5*0.01*12*(0.29690*sqrt(abs(p0(1)))-0.12600*p0(1)-0.35160*p0(1).^2+0.28430*p0(1).^3-0.10150*p0(1).^4)).^2;
    if abs(d1)<1e-4
        x = p0(1);
        y = 5*0.01*12*(0.29690*sqrt(abs(x))-0.12600*x-0.35160*x.^2+0.28430*x.^3-0.10150*x.^4);        
        p(i,2) = y*sign(p0(2));
    end
end

%parity = 0;
if parity==0
  t=[1,2,m+2; 1,m+2,m+1];
else
  t=[1,2,m+1; 2,m+2,m+1];
end

t=kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');
t=kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m);

% Reorder triangles in Cartesian order

ix=[];
for i=1:n-1
  ix1=i+(n-1)*(0:m-2);
  ix2=ix1+(n-1)*(m-1);

  if parity==0
    ix12=[ix2;ix1];
  else
    ix12=[ix1;ix2];
  end

  ix=[ix,reshape(ix12,1,[])];
end

t=t(ix,:);

