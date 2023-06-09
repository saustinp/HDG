function [f,f_udg] = sensor2d2(pg,udg,param,time)
%SENSOR2D2
%    [F,F_UDG] = SENSOR2D2(PG,UDG,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    29-Mar-2013 18:51:51
[ng,nc] = size(udg);
nch = 4;
nd = 2;
param8 = param{8};
param10 = param{10};
u1 = udg(:,1);
u2 = udg(:,2);
u3 = udg(:,3);
u5 = udg(:,5);
u6 = udg(:,6);
u9 = udg(:,9);
u11 = udg(:,11);
zero = zeros(ng,1);
t2 = 1.0./param8;
t3 = t2.*u1.*2.0e1;
t4 = exp(t3);
t5 = t4+2.0e1;
t6 = log(t5);
t7 = 1.0./t6;
t8 = 1.0./param10;
t9 = u6-t2.*t7.*u2.*u5.*2.0e1;
t10 = 1.0./param8.^2;
t11 = 1.0./t6.^2;
t12 = u11-t2.*t7.*u3.*u9.*2.0e1;
f = t8.*(t2.*t7.*t9.*2.0e1+t2.*t7.*t12.*2.0e1).*1.0e1-5.0;
if nargout > 1
    t13 = 1.0./t5;
    t14 = 1.0./param8.^3;
    t15 = 1.0./t6.^3;
    t16 = t2.*t7.*t8.*2.0e2;
    f_udg = [t8.*(t4.*t9.*t10.*t11.*t13.*4.0e2+t4.*t10.*t11.*t12.*t13.*4.0e2-t4.*t13.*t14.*t15.*u2.*u5.*8.0e3-t4.*t13.*t14.*t15.*u3.*u9.*8.0e3).*-1.0e1;t8.*t10.*t11.*u5.*-4.0e3;t8.*t10.*t11.*u9.*-4.0e3;zero;t8.*t10.*t11.*u2.*-4.0e3;t16;zero;zero;t8.*t10.*t11.*u3.*-4.0e3;zero;t16;zero];
end
% f = reshape(f,ng,nch,nd);
% f_udg = reshape(f_udg,ng,nch,nd,nc);
