function ui = ubouma5(p, nl, param, ib)

x = p(:,1);
y = p(:,2);
%ui = x.*nl(:,1) + y.*nl(:,2);
if ib==0
  ui = (x.^2 + y.^2)/2 - param(5);
elseif ib==1
  ui = x.*nl(:,1) + y.*nl(:,2);
end


