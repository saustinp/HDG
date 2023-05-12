function ui = ubou(p, nl, ug, param, ib)

%gb = fbou(xgf,nl,param,-ibf);
x = p(:,1);
y = p(:,2);
ui = sin(0.5*pi*x).*sin(0.5*pi*y);       
