function [sr,sr_udg] = source(p,udg,param,time)

[ng,nc] = size(udg);
nch = 4;

sr = zeros(ng,nch);
sr_udg = zeros(ng,nch,nc); 
