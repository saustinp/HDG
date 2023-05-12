function [sr, sr_udg] = sourcema3(p,udg,param,time)

x = p(:,1);
sr = -ones(size(x));

[ng,nc] = size(udg);
nch = 1;
sr_udg = zeros(ng,nch,nc); 
