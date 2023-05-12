function [f,f_udg] = avflux2d5(pg,udg,param,time)
%AVFLUX2D5
%    [F,F_UDG] = AVFLUX2D5(PG,UDG,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    29-Mar-2013 18:50:43
[ng,nc] = size(udg);
nch = 4;
nd = 2;
param1 = param{1};
param8 = param{8};
param10 = param{10};
u1 = udg(:,1);
u2 = udg(:,2);
u3 = udg(:,3);
u5 = udg(:,5);
u6 = udg(:,6);
u7 = udg(:,7);
u8 = udg(:,8);
u9 = udg(:,9);
u10 = udg(:,10);
u11 = udg(:,11);
u12 = udg(:,12);
zero = zeros(ng,1);
t2 = 1.0./param8;
t3 = t2.*u1.*2.0e1;
t4 = exp(t3);
t5 = t4+2.0e1;
t6 = log(t5);
t7 = 1.0./t6;
t8 = 1.0./param10;
t16 = t2.*t7.*u2.*u5.*2.0e1;
t9 = -t16+u6;
t10 = t2.*t7.*t9.*2.0e1;
t17 = t2.*t7.*u3.*u9.*2.0e1;
t11 = -t17+u11;
t12 = t2.*t7.*t11.*2.0e1;
t13 = t10+t12;
t14 = t8.*t13;
t15 = t14-1.0./2.0;
t18 = 1.0./param8.^2;
t19 = 1.0./t6.^2;
t20 = param1-1.0;
t21 = 1.0./t5;
t22 = u2.^2;
t23 = t18.*t19.*t22.*2.0e2;
t24 = u3.^2;
t25 = t18.*t19.*t24.*2.0e2;
t26 = t23+t25;
t27 = 1.0./param8.^3;
t28 = 1.0./t6.^3;
t29 = t4.*t9.*t18.*t19.*t21.*4.0e2;
t30 = t4.*t11.*t18.*t19.*t21.*4.0e2;
t32 = t4.*t21.*t27.*t28.*u2.*u5.*8.0e3;
t33 = t4.*t21.*t27.*t28.*u3.*u9.*8.0e3;
t31 = t29+t30-t32-t33;
t34 = t9.*t18.*t19.*u2.*4.0e2;
t42 = t2.*t7.*u3.*u5.*2.0e1;
t35 = -t42+u7;
t36 = t18.*t19.*t35.*u3.*4.0e2;
t37 = t34+t36;
t38 = param8.*t6.*t37.*(1.0./2.0e1);
t39 = t4.*t21.*t26.*u5;
t40 = t38+t39-u8;
t60 = t20.*t40;
t41 = -t60+u8;
t43 = 1.0./param8.^4;
t44 = 1.0./t6.^4;
t53 = t2.*t7.*u2.*u9.*2.0e1;
t45 = -t53+u10;
t46 = t18.*t19.*t45.*u2.*4.0e2;
t47 = t11.*t18.*t19.*u3.*4.0e2;
t48 = t46+t47;
t49 = param8.*t6.*t48.*(1.0./2.0e1);
t50 = t4.*t21.*t26.*u9;
t51 = t49+t50-u12;
t61 = t20.*t51;
t52 = -t61+u12;
f = [t15.*u5;t15.*u6;t15.*u7;t15.*t41;t15.*u9;t15.*u10;t15.*u11;t15.*t52];
if nargout > 1
    t54 = t4.*t21.*t22.*t27.*t28.*8.0e3;
    t55 = t4.*t21.*t24.*t27.*t28.*8.0e3;
    t56 = t54+t55;
    t57 = t2.*u1.*4.0e1;
    t58 = exp(t57);
    t59 = 1.0./t5.^2;
    t62 = t22.*t27.*t28.*8.0e3;
    t63 = t24.*t27.*t28.*8.0e3;
    t64 = t62+t63;
    t65 = param8.*t6.*t64.*(1.0./2.0e1);
    t66 = t15.*t20.*(t65-t4.*t21.*t26);
    t67 = t2.*t7.*t8.*u5.*2.0e1;
    t68 = t2.*t7.*t8.*u6.*2.0e1;
    t69 = t2.*t7.*t8.*u7.*2.0e1;
    t70 = t2.*t7.*t8.*t41.*2.0e1;
    t71 = t2.*t7.*t8.*u9.*2.0e1;
    t72 = t2.*t7.*t8.*u10.*2.0e1;
    t73 = t2.*t7.*t8.*u11.*2.0e1;
    t74 = t2.*t7.*t8.*t52.*2.0e1;
    t75 = param1.*t15;
    f_udg = [-t8.*t31.*u5;-t8.*t31.*u6;-t8.*t31.*u7;-t8.*t31.*t41+t15.*t20.*(param8.*t6.*(t4.*t9.*t21.*t27.*t28.*u2.*1.6e4+t4.*t21.*t27.*t28.*t35.*u3.*1.6e4-t4.*t21.*t22.*t43.*t44.*u5.*1.6e5-t4.*t21.*t24.*t43.*t44.*u5.*1.6e5).*(1.0./2.0e1)-t4.*t21.*t37+t4.*t21.*t56.*u5-t2.*t4.*t21.*t26.*u5.*2.0e1+t2.*t26.*t58.*t59.*u5.*2.0e1);-t8.*t31.*u9;-t8.*t31.*u10;-t8.*t31.*u11;-t8.*t31.*t52+t15.*t20.*(param8.*t6.*(t4.*t11.*t21.*t27.*t28.*u3.*1.6e4+t4.*t21.*t27.*t28.*t45.*u2.*1.6e4-t4.*t21.*t22.*t43.*t44.*u9.*1.6e5-t4.*t21.*t24.*t43.*t44.*u9.*1.6e5).*(1.0./2.0e1)-t4.*t21.*t48+t4.*t21.*t56.*u9-t2.*t4.*t21.*t26.*u9.*2.0e1+t2.*t26.*t58.*t59.*u9.*2.0e1);t8.*t18.*t19.*u5.^2.*-4.0e2;t8.*t18.*t19.*u5.*u6.*-4.0e2;t8.*t18.*t19.*u5.*u7.*-4.0e2;-t15.*t20.*(param8.*t6.*(t9.*t18.*t19.*4.0e2-t27.*t28.*u2.*u5.*8.0e3).*(1.0./2.0e1)+t4.*t18.*t19.*t21.*u2.*u5.*4.0e2)-t8.*t18.*t19.*t41.*u5.*4.0e2;t8.*t18.*t19.*u5.*u9.*-4.0e2;t8.*t18.*t19.*u5.*u10.*-4.0e2;t8.*t18.*t19.*u5.*u11.*-4.0e2;-t15.*t20.*(param8.*t6.*(t18.*t19.*t45.*4.0e2-t27.*t28.*u2.*u9.*8.0e3).*(1.0./2.0e1)+t4.*t18.*t19.*t21.*u2.*u9.*4.0e2)-t8.*t18.*t19.*t52.*u5.*4.0e2;t8.*t18.*t19.*u5.*u9.*-4.0e2;t8.*t18.*t19.*u6.*u9.*-4.0e2;t8.*t18.*t19.*u7.*u9.*-4.0e2;-t15.*t20.*(param8.*t6.*(t18.*t19.*t35.*4.0e2-t27.*t28.*u3.*u5.*8.0e3).*(1.0./2.0e1)+t4.*t18.*t19.*t21.*u3.*u5.*4.0e2)-t8.*t18.*t19.*t41.*u9.*4.0e2;t8.*t18.*t19.*u9.^2.*-4.0e2;t8.*t18.*t19.*u9.*u10.*-4.0e2;t8.*t18.*t19.*u9.*u11.*-4.0e2;-t15.*t20.*(param8.*t6.*(t11.*t18.*t19.*4.0e2-t27.*t28.*u3.*u9.*8.0e3).*(1.0./2.0e1)+t4.*t18.*t19.*t21.*u3.*u9.*4.0e2)-t8.*t18.*t19.*t52.*u9.*4.0e2;zero;zero;zero;zero;zero;zero;zero;zero;t14-t8.*t18.*t19.*u2.*u5.*4.0e2-1.0./2.0;t8.*t18.*t19.*u2.*u6.*-4.0e2;t8.*t18.*t19.*u2.*u7.*-4.0e2;t66-t8.*t18.*t19.*t41.*u2.*4.0e2;t8.*t18.*t19.*u2.*u9.*-4.0e2;t8.*t18.*t19.*u2.*u10.*-4.0e2;t8.*t18.*t19.*u2.*u11.*-4.0e2;t8.*t18.*t19.*t52.*u2.*-4.0e2;t67;t14+t68-1.0./2.0;t69;t70-t2.*t7.*t15.*t20.*u2.*2.0e1;t71;t72;t73;t74;zero;zero;t15;t2.*t7.*t15.*t20.*u3.*-2.0e1;zero;zero;zero;zero;zero;zero;zero;t75;zero;zero;zero;zero;t8.*t18.*t19.*u3.*u5.*-4.0e2;t8.*t18.*t19.*u3.*u6.*-4.0e2;t8.*t18.*t19.*u3.*u7.*-4.0e2;t8.*t18.*t19.*t41.*u3.*-4.0e2;t14-t8.*t18.*t19.*u3.*u9.*4.0e2-1.0./2.0;t8.*t18.*t19.*u3.*u10.*-4.0e2;t8.*t18.*t19.*u3.*u11.*-4.0e2;t66-t8.*t18.*t19.*t52.*u3.*4.0e2;zero;zero;zero;zero;zero;t15;zero;t2.*t7.*t15.*t20.*u2.*-2.0e1;t67;t68;t69;t70;t71;t72;t14+t73-1.0./2.0;t74-t2.*t7.*t15.*t20.*u3.*2.0e1;zero;zero;zero;zero;zero;zero;zero;t75];
end
f = reshape(f,ng,nch,nd);
f_udg = reshape(f_udg,ng,nch,nd,nc);
