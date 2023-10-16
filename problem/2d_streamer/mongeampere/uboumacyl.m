function [ui, dui] = uboumacyl(p, nl, param, ib, udg)

x = p(:,1);
y = p(:,2);

dui = 0*p;

if ib==0
  ui = (x.^2 + y.^2)/2 - param(5);
elseif ib==1
  ui = x.*nl(:,1) + y.*nl(:,2);
elseif ib==2
  ui = -ones(size(x));  
end

if (ib==2) && nargin>4
  ui = 1 - udg(:,1).*udg(:,1)/param(1) - udg(:,2).*udg(:,2)/param(2);      
  dui(:,1) = -2*udg(:,1)/param(1);
  dui(:,2) = -2*udg(:,2)/param(2);      
end
