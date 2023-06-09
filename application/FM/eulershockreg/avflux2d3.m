function [f,f_udg] = avflux2d3(pg,udg,param,time)
%AVFLUX2D3
%    [F,F_UDG] = AVFLUX2D3(PG,UDG,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    29-Mar-2013 18:50:18
[ng,nc] = size(udg);
nch = 4;
nd = 2;
param1 = param{1};
param9 = param{9};
param10 = param{10};
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
zero = zeros(ng,1);
t2 = 1.0./u1;
t3 = 1.0./param10;
t11 = t2.*u2.*u5;
t4 = -t11+u6;
t5 = t2.*t4;
t12 = t2.*u3.*u9;
t6 = -t12+u11;
t7 = t2.*t6;
t8 = t5+t7;
t9 = t3.*t8;
t10 = t9-1.0./2.0;
t13 = 1.0./u1.^2;
t14 = param1-1.0;
t15 = 1.0./param9;
t16 = u2.^2;
t17 = t13.*t16.*(1.0./2.0);
t18 = u3.^2;
t19 = t13.*t18.*(1.0./2.0);
t20 = t17+t19;
t24 = t20.*u1;
t21 = -t24+u4;
t22 = t14.*t15.*t21.*2.0e1;
t23 = exp(t22);
t25 = t23+2.0e1;
t26 = 1.0./t25;
t27 = 1.0./u1.^3;
t28 = t4.*t13;
t29 = t6.*t13;
t31 = t27.*u2.*u5;
t32 = t27.*u3.*u9;
t30 = t28+t29-t31-t32;
t33 = 1.0./u1.^4;
t36 = t2.*u3.*u5;
t34 = -t36+u7;
t35 = t4.*t13.*u2;
t37 = t13.*t34.*u3;
t38 = t16.*t27;
t39 = t18.*t27;
t40 = t38+t39;
t41 = t20.*u5;
t42 = t35+t37;
t43 = t42.*u1;
t44 = t41+t43-u8;
t45 = t14.^2;
t52 = t40.*u1;
t46 = t17+t19-t52;
t61 = t14.*t23.*t26.*t44;
t47 = -t61+u8;
t49 = t2.*u2.*u9;
t48 = -t49+u10;
t50 = t13.*t48.*u2;
t51 = t6.*t13.*u3;
t53 = t20.*u9;
t54 = t50+t51;
t55 = t54.*u1;
t56 = t53+t55-u12;
t57 = t14.*t15.*t21.*4.0e1;
t58 = exp(t57);
t59 = 1.0./t25.^2;
t62 = t14.*t23.*t26.*t56;
t60 = -t62+u12;
f = [t10.*u5;t10.*u6;t10.*u7;t10.*t47;t10.*u9;t10.*u10;t10.*u11;t10.*t60];
if nargout > 1
    t63 = t2.*t3.*u5;
    t64 = t2.*t3.*u6;
    t65 = t2.*t3.*u7;
    t66 = t2.*t3.*t47;
    t67 = t2.*t3.*u9;
    t68 = t2.*t3.*u10;
    t69 = t2.*t3.*u11;
    t70 = t2.*t3.*t60;
    t71 = t14.*t23.*t26;
    t72 = t71+1.0;
    t73 = t10.*t72;
    f_udg = [-t3.*t30.*u5;-t3.*t30.*u6;-t3.*t30.*u7;-t10.*(t14.*t23.*t26.*(t35+t37-u1.*(t4.*t27.*u2.*2.0-t16.*t33.*u5-t18.*t33.*u5+t27.*t34.*u3.*2.0)-t40.*u5)-t15.*t23.*t26.*t44.*t45.*t46.*2.0e1+t15.*t44.*t45.*t46.*t58.*t59.*2.0e1)-t3.*t30.*t47;-t3.*t30.*u9;-t3.*t30.*u10;-t3.*t30.*u11;-t10.*(t14.*t23.*t26.*(t50+t51-u1.*(t6.*t27.*u3.*2.0-t16.*t33.*u9-t18.*t33.*u9+t27.*t48.*u2.*2.0)-t40.*u9)-t15.*t23.*t26.*t45.*t46.*t56.*2.0e1+t15.*t45.*t46.*t56.*t58.*t59.*2.0e1)-t3.*t30.*t60;-t3.*t13.*u5.^2;-t3.*t13.*u5.*u6;-t3.*t13.*u5.*u7;-t10.*(t14.*t23.*t26.*(u1.*(t28-t31)+t13.*u2.*u5)-t2.*t15.*t23.*t26.*t44.*t45.*u2.*2.0e1+t2.*t15.*t44.*t45.*t58.*t59.*u2.*2.0e1)-t3.*t13.*t47.*u5;-t3.*t13.*u5.*u9;-t3.*t13.*u5.*u10;-t3.*t13.*u5.*u11;-t10.*(t14.*t23.*t26.*(u1.*(t13.*t48-t27.*u2.*u9)+t13.*u2.*u9)-t2.*t15.*t23.*t26.*t45.*t56.*u2.*2.0e1+t2.*t15.*t45.*t56.*t58.*t59.*u2.*2.0e1)-t3.*t13.*t60.*u5;-t3.*t13.*u5.*u9;-t3.*t13.*u6.*u9;-t3.*t13.*u7.*u9;-t10.*(t14.*t23.*t26.*(u1.*(t13.*t34-t27.*u3.*u5)+t13.*u3.*u5)-t2.*t15.*t23.*t26.*t44.*t45.*u3.*2.0e1+t2.*t15.*t44.*t45.*t58.*t59.*u3.*2.0e1)-t3.*t13.*t47.*u9;-t3.*t13.*u9.^2;-t3.*t13.*u9.*u10;-t3.*t13.*u9.*u11;-t10.*(t14.*t23.*t26.*(u1.*(t29-t32)+t13.*u3.*u9)-t2.*t15.*t23.*t26.*t45.*t56.*u3.*2.0e1+t2.*t15.*t45.*t56.*t58.*t59.*u3.*2.0e1)-t3.*t13.*t60.*u9;zero;zero;zero;-t10.*(t15.*t23.*t26.*t44.*t45.*2.0e1-t15.*t44.*t45.*t58.*t59.*2.0e1);zero;zero;zero;-t10.*(t15.*t23.*t26.*t45.*t56.*2.0e1-t15.*t45.*t56.*t58.*t59.*2.0e1);t9-t3.*t13.*u2.*u5-1.0./2.0;-t3.*t13.*u2.*u6;-t3.*t13.*u2.*u7;-t3.*t13.*t47.*u2-t10.*t14.*t23.*t26.*t46;-t3.*t13.*u2.*u9;-t3.*t13.*u2.*u10;-t3.*t13.*u2.*u11;-t3.*t13.*t60.*u2;t63;t9+t64-1.0./2.0;t65;t66-t2.*t10.*t14.*t23.*t26.*u2;t67;t68;t69;t70;zero;zero;t10;-t2.*t10.*t14.*t23.*t26.*u3;zero;zero;zero;zero;zero;zero;zero;t73;zero;zero;zero;zero;-t3.*t13.*u3.*u5;-t3.*t13.*u3.*u6;-t3.*t13.*u3.*u7;-t3.*t13.*t47.*u3;t9-t3.*t13.*u3.*u9-1.0./2.0;-t3.*t13.*u3.*u10;-t3.*t13.*u3.*u11;-t3.*t13.*t60.*u3-t10.*t14.*t23.*t26.*t46;zero;zero;zero;zero;zero;t10;zero;-t2.*t10.*t14.*t23.*t26.*u2;t63;t64;t65;t66;t67;t68;t9+t69-1.0./2.0;t70-t2.*t10.*t14.*t23.*t26.*u3;zero;zero;zero;zero;zero;zero;zero;t73];
end
f = reshape(f,ng,nch,nd);
f_udg = reshape(f_udg,ng,nch,nd,nc);
