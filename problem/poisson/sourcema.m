function sr = sourcema(p,udg,param,time)

x = p(:,1);
y = p(:,2);
sr = exp(x.^2 + y.^2).*(x.^2 + y.^2 + 1);
%sr = -exp(x.^2/2 + y.^2/2).*(x.^2 + y.^2 + 2);


