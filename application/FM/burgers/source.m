
function [sr,sr_udg] = source(p,udg,param,time)

[ng,nch] = size(udg);
nc = nch;
nd = nch;

if nd ~= 1; error('Routine not valid for nd > 1.'); end
warning('Is cos() correct in the source term?');

global Af;
global Nc;    % Cut-off
global sigma_N;
global deltaT;
global L;
% global ne;
% global porder;
% global Phi_jXi;
% global alpha_j;
% global beta_j;
% global gpvl;
global hackToSpeed;

% h = L/ne;
% xc = - L/2 + h/2 + h*(0:(ne-1))';

sr = zeros(ng,1);

implementationFlag = 1;
if implementationFlag == 0      % Inexact integration
    for i=(-Nc:Nc)
        sr = sr + sigma_N(i+Nc+1) * cos(2*pi*i*p(:) / L) / sqrt(max(abs(i),1));
    end
    sr = (Af / sqrt(deltaT)) * sr;
    
elseif implementationFlag == 1 % Exact integration of Galerkin projection
    for i=(-Nc:Nc)
        sr = sr + (sigma_N(i+Nc+1) / sqrt(max(abs(i),1))) .* real(hackToSpeed{i+Nc+1});
    end
    sr = (Af / sqrt(deltaT)) * sr;
    
elseif implementationFlag == 2  %Exact integration of Galerkin + Fourier projection
    error('Exact integration of Galerkin + Fourier projection not implemented yet.');
    % Just replace alpha_j -> beta_j
end

% % % sr = - (2*pi/L)*cos(2*pi*p(:)/L).*sin(2*pi*p(:)/L);

sr_udg = zeros(ng,nc,nc);

if any(isnan(sr(:))) || any(isnan(sr_udg(:))); error('NaN detected'); end

end
