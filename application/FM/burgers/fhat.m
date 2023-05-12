
function [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time)

alpha = 1.0e4; % Parameter for smooth surrogate of abs() function
beta = param{1};    % Upwind parameter

[ng,nch] = size(uh);
nd = nch;

if nd ~= 1; error('Routine not valid for nd > 1.'); end

u = udg(:,1:nch);

[f,f_uhdg] = flux(p,uh,param,time);

tauMin = 1.0e-6;%0.1;

An1  = zeros(ng,1,1);
An2  = zeros(ng,1,1);
An1_uh = zeros(ng,1,1,1);
An2_u = zeros(ng,1,1,1);
tmp1 = uh(:,1) .* nl(:,1);
tmp2 = u(:,1) .* nl(:,1);
An1(:,1,1) = beta * ( 2*(tmp1.*(atan(alpha*tmp1)/pi + 1/2) + 1/2 - atan(alpha)/pi) - tmp1 );    % An(:,1,1) = beta * abs(tmp);
An1_uh(:,1,1,1) = beta * nl(:,1) .* ( 2*(atan(alpha*tmp1)/pi + (alpha*tmp1)./(pi*(alpha^2*tmp1.^2 + 1)) + 1/2) - 1 );   % An_uh(:,1,1,1) = beta * nl(:,1) .* sign(tmp);
An2(:,1,1) = 1*beta * ( 2*(tmp2.*(atan(alpha*tmp2)/pi + 1/2) + 1/2 - atan(alpha)/pi) - tmp2 );
An2_u(:,1,1,1) = 1*beta * nl(:,1) .* ( 2*(atan(alpha*tmp2)/pi + (alpha*tmp2)./(pi*(alpha^2*tmp2.^2 + 1)) + 1/2) - 1 );

An = tauMin + An1 + An2;
An_uh = An1_uh;
An_u = An2_u;

fh = f .* nl + An .* (u-uh);

fh_udg = An_u .* (u-uh) + An;
fh_uh = f_uhdg .* nl + An_uh .* (u-uh) - An;

if any(isnan(fh(:))) || any(isnan(fh_udg(:))) || any(isnan(fh_uh(:))); error('NaN detected'); end

end

function [y,dydx] = smoothAbs(x) 

alpha = 10;
y = 2*(x.*(atan(alpha*x)/pi + 1/2) + 1/2 - atan(alpha)/pi) - x;
dydx = 2*(atan(alpha*x)/pi + (alpha*x)./(pi*(alpha^2*x.^2 + 1)) + 1/2) - 1;

end
