function [f,f_udg] = sensor2d1(pg,udg,param,time)
%SENSOR2D1
%    [F,F_UDG] = SENSOR2D1(PG,UDG,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    31-Mar-2013 00:15:36
[ng,nc] = size(udg);
nch = 4;
nd = 2;
param10 = param{10};
u1 = udg(:,1);
u2 = udg(:,2);
u3 = udg(:,3);
u5 = udg(:,5);
u6 = udg(:,6);
u9 = udg(:,9);
u11 = udg(:,11);
zero = zeros(ng,1);
t2 = 1.0./u1;
t3 = 1.0./param10;
t4 = u6-t2.*u2.*u5;
t5 = 1.0./u1.^2;
t6 = u11-t2.*u3.*u9;
f = t3.*(t2.*t4+t2.*t6).*1.0e1-5.0;
if nargout > 1
    t7 = 1.0./u1.^3;
    t8 = t2.*t3.*1.0e1;
    f_udg = [t3.*(t4.*t5+t5.*t6-t7.*u2.*u5-t7.*u3.*u9).*-1.0e1;t3.*t5.*u5.*-1.0e1;t3.*t5.*u9.*-1.0e1;zero;t3.*t5.*u2.*-1.0e1;t8;zero;zero;t3.*t5.*u3.*-1.0e1;zero;t8;zero];
end
% f = reshape(f,ng,nch,nd);
% f_udg = reshape(f_udg,ng,nch,nd,nc);
