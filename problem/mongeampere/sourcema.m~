function [sr, sr_udg] = sourcema(p,udg,param,time)

x = p(:,1);
y = p(:,2);
sr = -exp(x.^2 + y.^2).*(x.^2 + y.^2 + 1);
%sr = -exp(x.^2/2 + y.^2/2).*(x.^2 + y.^2 + 2);

% if param(end) == 2
%   uex = exp(0.5*(x.*x+y.*y));
%   uxx = uex + x.^2.*uex;
%   uyy = uex + y.^2.*uex;  
%   sr = -(uxx+uyy);
% elseif param(end) == 3
%   sr = -(udg(:,2) + udg(:,3));
% end
  
[ng,nc] = size(udg);
nch = 1;
sr_udg = zeros(ng,nch,nc); 

