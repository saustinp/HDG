function ui = ubouma(p, nl, param, ib)

x = p(:,1);
y = p(:,2);
ui = exp(0.5*(x.*x+y.*y));       
