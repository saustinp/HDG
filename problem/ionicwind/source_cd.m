function [sr,sr_udg] = source_cd(p,udg,param,time)

[ng,nc] = size(udg);
nch = 1;

sr = zeros(ng,nch); 
sr_udg = zeros(ng,nch,nc); 



