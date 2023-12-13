function [f,f_udg] = flux2d(p,udg,param1,time)
%FLUX2D
%    [F,F_UDG] = FLUX2D(P,UDG,PARAM1,TIME)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    06-Sep-2023 11:52:06
[ng,nc] = size(udg);
nch = 1;
nd = 2;
one = ones(ng,1);
u1 = udg(:,1);
u2 = udg(:,2);
u3 = udg(:,3);
x1 = p(:,1);
zero = zeros(ng,1);
f = [(u2.*x1)./1.0e+16;x1.*(u1+u3./1.0e+16)];
if nargout > 1
    t2 = (one.*x1)./1.0e+16;
    f_udg = [zero;one.*x1;t2;zero;zero;t2];
end

f = reshape(f,ng,nch,nd);
f_udg = reshape(f_udg,ng,nch,nd,nc);