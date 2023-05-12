function [fh,fh_udg,fh_uh] = fhat_ns2d(nl,pg,udg,uh,param,time)
%FHAT_NS2D
%    [FH,FH_UDG,FH_UH] = FHAT_NS2D(NL,PG,UDG,UH,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    21-Dec-2019 10:43:32
[ng,nc] = size(udg);
nch = 4;
nd = 2;
nl1 = nl(:,1);
nl2 = nl(:,2);
one = ones(ng,1);
param1 = param{1};
param2 = param{2};
param3 = param{3};
param4 = param{4};
param5 = param{5};
param6 = param{6};
u1 = udg(:,1);
u2 = udg(:,2);
u3 = udg(:,3);
u4 = udg(:,4);
u5 = udg(:,5);
u6 = udg(:,6);
u7 = udg(:,7);
u8 = udg(:,8);
u9 = udg(:,9);
u10 = udg(:,10);
u11 = udg(:,11);
u12 = udg(:,12);
uh4 = uh(:,4);
x1 = pg(:,1);
x2 = pg(:,2);
zero = zeros(ng,1);
t2 = param2.^2;
t3 = 1.0./u1.^2;
t4 = 1.0./u1;
t5 = param1-1.0;
t6 = x1.^2;
t7 = t2.*t6.*(1.0./2.0);
t8 = x2.^2;
t9 = t2.*t8.*(1.0./2.0);
t10 = t7+t9;
t11 = u2.^2;
t12 = t3.*t11.*(1.0./2.0);
t13 = u3.^2;
t14 = t3.*t13.*(1.0./2.0);
t15 = t12+t14;
t16 = param5.^2;
t17 = t10.*u1;
t33 = t15.*u1;
t18 = t17-t33+u4;
t19 = 1.0./param3;
t20 = 1.0./param4;
t21 = 1.0./param5.^2;
t22 = 1.0./t5;
t23 = t15.*u5;
t50 = t4.*u2.*u5;
t24 = -t50+u6;
t25 = t3.*t24.*u2;
t51 = t4.*u3.*u5;
t26 = -t51+u7;
t27 = t3.*t26.*u3;
t28 = t25+t27;
t29 = t28.*u1;
t30 = t2.*u1.*x1;
t49 = t10.*u5;
t31 = t23+t29+t30-t49-u8;
t32 = t5.*t31.*u1;
t34 = t5.*t18.*u5;
t35 = t32+t34;
t36 = 1.0./u1.^3;
t37 = t15.*u9;
t57 = t4.*u2.*u9;
t38 = -t57+u10;
t39 = t3.*t38.*u2;
t58 = t4.*u3.*u9;
t40 = -t58+u11;
t41 = t3.*t40.*u3;
t42 = t39+t41;
t43 = t42.*u1;
t44 = t2.*u1.*x2;
t56 = t10.*u9;
t45 = t37+t43+t44-t56-u12;
t46 = t5.*t45.*u1;
t47 = t5.*t18.*u9;
t48 = t46+t47;
fh = param6.*(u4-uh4)-t19.*t20.*t21.*t22.*(nl1.*param1.*t3.*t16.*t35+nl2.*param1.*t3.*t16.*t48);
if nargout > 1
    t52 = t11.*t36;
    t53 = t13.*t36;
    t54 = t52+t53;
    t55 = 1.0./u1.^4;
    t59 = t54.*u1;
    t60 = t7+t9-t12-t14+t59;
    t61 = t5.*t18;
    t62 = t61-t5.*t60.*u1;
    fh_udg = [t19.*t20.*t21.*t22.*(nl1.*param1.*t16.*t35.*t36.*2.0+nl2.*param1.*t16.*t36.*t48.*2.0-nl1.*param1.*t3.*t16.*(t5.*t31+t5.*u1.*(t25+t27-u1.*(t24.*t36.*u2.*2.0+t26.*t36.*u3.*2.0-t11.*t55.*u5-t13.*t55.*u5)-t54.*u5+t2.*x1)+t5.*t60.*u5)-nl2.*param1.*t3.*t16.*(t5.*t45+t5.*u1.*(t39+t41+u1.*(t11.*t55.*u9-t36.*t38.*u2.*2.0+t13.*t55.*u9-t36.*t40.*u3.*2.0)-t54.*u9+t2.*x2)+t5.*t60.*u9));-t19.*t20.*t21.*t22.*(nl1.*param1.*t3.*t16.*(t5.*u1.*(u1.*(t3.*t24-t36.*u2.*u5)+t3.*u2.*u5)-t4.*t5.*u2.*u5)+nl2.*param1.*t3.*t16.*(t5.*u1.*(u1.*(t3.*t38-t36.*u2.*u9)+t3.*u2.*u9)-t4.*t5.*u2.*u9));-t19.*t20.*t21.*t22.*(nl1.*param1.*t3.*t16.*(t5.*u1.*(u1.*(t3.*t26-t36.*u3.*u5)+t3.*u3.*u5)-t4.*t5.*u3.*u5)+nl2.*param1.*t3.*t16.*(t5.*u1.*(u1.*(t3.*t40-t36.*u3.*u9)+t3.*u3.*u9)-t4.*t5.*u3.*u9));param6-t19.*t20.*t21.*t22.*(nl1.*param1.*t3.*t5.*t16.*u5+nl2.*param1.*t3.*t5.*t16.*u9);-nl1.*param1.*t3.*t19.*t20.*t22.*t62;-nl1.*param1.*t3.*t19.*t20.*u2;-nl1.*param1.*t3.*t19.*t20.*u3;nl1.*param1.*t4.*t19.*t20;-nl2.*param1.*t3.*t19.*t20.*t22.*t62;-nl2.*param1.*t3.*t19.*t20.*u2;-nl2.*param1.*t3.*t19.*t20.*u3;nl2.*param1.*t4.*t19.*t20];
end
if nargout > 2
    fh_uh = [zero;zero;zero;-one.*param6];
end
fh = reshape(fh,ng,1);
fh_udg = reshape(fh_udg,ng,1,nc);
fh_uh = reshape(fh_uh,ng,1,nch);
