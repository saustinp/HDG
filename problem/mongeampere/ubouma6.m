function ui = ubouma6(p, nl, param, ib)

x = p(:,1);
y = p(:,2);
%ui = x.*nl(:,1) + y.*nl(:,2);
if ib==0
  ui = (x.^2 + y.^2)/2 - param(5);
elseif ib==1
  ui = x.*nl(:,1) + y.*nl(:,2);
elseif ib==2
  ui = x.*x + y.*y;  
%   [p nl]
%   ui = 0.5*((x.*x + nl(:,1).*nl(:,1)) + (y.*y + nl(:,2).*nl(:,2)));  
%   pause
end


