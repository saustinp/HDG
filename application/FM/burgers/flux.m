
function [f,f_udg] = flux(p,udg,param,time)

[ng,nch] = size(udg);
nd = nch;

if nd ~= 1; error('Routine not valid for nd > 1.'); end

f = zeros(ng,1,1);
f_udg = zeros(ng,1,1,1);
f(:,1,1) = udg(:,1).^2 / 2;
f_udg(:,1,1,1) = udg(:,1);

if any(isnan(f(:))) || any(isnan(f_udg(:))); error('NaN detected'); end

end
