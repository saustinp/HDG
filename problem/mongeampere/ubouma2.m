function ui = ubouma2(p, nl, param, ib)

R = param(1);

x = p(:,1);
y = p(:,2);
ui = -sqrt(R^2 - x.^2 - y.^2);  
