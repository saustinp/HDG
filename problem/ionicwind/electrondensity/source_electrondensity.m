function [sr,sr_udg] = source_electrondensity(p,udg,param,time)

[ng,nc] = size(udg);
nch = 1;

% Physics parameters
% r0 = param{1};
% z0 = param{2};
% s0 = param{3};
% Nmax = param{4};
% e = param{5};
% epsilon0 = param{6};
% Ua = param{7};
% gamma = param{8};
% E_bd = param{9};
% r_tip = param{10};
% n_ref = param{11};
% N = param{12};
% mue_ref = param{13};
% D_star = param{14};

sr = zeros(ng,nch); 
sr_udg = zeros(ng,nch,nc);