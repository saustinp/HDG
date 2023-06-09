function [f,f_udg] = avflux2d6(pg,udg,param,time)
%AVFLUX2D6
%    [F,F_UDG] = AVFLUX2D6(PG,UDG,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    29-Mar-2013 18:51:01
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
t19 = t2.*t7.*u2.*u5.*2.0e1;
t9 = -t19+u6;
t10 = t2.*t7.*t9.*2.0e1;
t20 = t2.*t7.*u3.*u9.*2.0e1;
t11 = -t20+u11;
t12 = t2.*t7.*t11.*2.0e1;
t13 = t10+t12;
t14 = t8.*t13.*1.0e1;
t15 = t14-5.0;
t16 = exp(t15);
t17 = t16+1.0;
t18 = log(t17);
t21 = 1.0./param8.^2;
t22 = 1.0./t6.^2;
t23 = param1-1.0;
t24 = 1.0./t5;
t25 = u2.^2;
t26 = t21.*t22.*t25.*2.0e2;
t27 = u3.^2;
t28 = t21.*t22.*t27.*2.0e2;
t29 = t26+t28;
t30 = 1.0./param8.^3;
t31 = 1.0./t6.^3;
t32 = 1.0./t17;
t33 = t4.*t9.*t21.*t22.*t24.*4.0e2;
t34 = t4.*t11.*t21.*t22.*t24.*4.0e2;
t36 = t4.*t24.*t30.*t31.*u2.*u5.*8.0e3;
t37 = t4.*t24.*t30.*t31.*u3.*u9.*8.0e3;
t35 = t33+t34-t36-t37;
t38 = t9.*t21.*t22.*u2.*4.0e2;
t44 = t2.*t7.*u3.*u5.*2.0e1;
t39 = -t44+u7;
t40 = t21.*t22.*t39.*u3.*4.0e2;
t41 = t38+t40;
t42 = 1.0./param8.^4;
t43 = 1.0./t6.^4;
t45 = param8.*t6.*t41.*(1.0./2.0e1);
t46 = t4.*t24.*t29.*u5;
t47 = t45+t46-u8;
t64 = t23.*t47;
t48 = -t64+u8;
t53 = t2.*t7.*u2.*u9.*2.0e1;
t49 = -t53+u10;
t50 = t21.*t22.*t49.*u2.*4.0e2;
t51 = t11.*t21.*t22.*u3.*4.0e2;
t52 = t50+t51;
t54 = t4.*t24.*t25.*t30.*t31.*8.0e3;
t55 = t4.*t24.*t27.*t30.*t31.*8.0e3;
t56 = t54+t55;
t57 = t2.*u1.*4.0e1;
t58 = exp(t57);
t59 = 1.0./t5.^2;
t60 = param8.*t6.*t52.*(1.0./2.0e1);
t61 = t4.*t24.*t29.*u9;
t62 = t60+t61-u12;
t65 = t23.*t62;
t63 = -t65+u12;
f = [t18.*u5.*(1.0./1.0e1);t18.*u6.*(1.0./1.0e1);t18.*u7.*(1.0./1.0e1);t18.*t48.*(1.0./1.0e1);t18.*u9.*(1.0./1.0e1);t18.*u10.*(1.0./1.0e1);t18.*u11.*(1.0./1.0e1);t18.*t63.*(1.0./1.0e1)];
if nargout > 1
    t66 = t18.*(1.0./1.0e1);
    t67 = t25.*t30.*t31.*8.0e3;
    t68 = t27.*t30.*t31.*8.0e3;
    t69 = t67+t68;
    t70 = param8.*t6.*t69.*(1.0./2.0e1);
    t71 = t18.*t23.*(t70-t4.*t24.*t29).*(1.0./1.0e1);
    t72 = t2.*t7.*t8.*t16.*t32.*u5.*2.0e1;
    t73 = t2.*t7.*t8.*t16.*t32.*u6.*2.0e1;
    t74 = t2.*t7.*t8.*t16.*t32.*u7.*2.0e1;
    t75 = t2.*t7.*t8.*t16.*t32.*t48.*2.0e1;
    t76 = t2.*t7.*t8.*t16.*t32.*u9.*2.0e1;
    t77 = t2.*t7.*t8.*t16.*t32.*u10.*2.0e1;
    t78 = t2.*t7.*t8.*t16.*t32.*u11.*2.0e1;
    t79 = t2.*t7.*t8.*t16.*t32.*t63.*2.0e1;
    t80 = param1.*t18.*(1.0./1.0e1);
    f_udg = [-t8.*t16.*t32.*t35.*u5;-t8.*t16.*t32.*t35.*u6;-t8.*t16.*t32.*t35.*u7;t18.*t23.*(param8.*t6.*(t4.*t9.*t24.*t30.*t31.*u2.*1.6e4+t4.*t24.*t30.*t31.*t39.*u3.*1.6e4-t4.*t24.*t25.*t42.*t43.*u5.*1.6e5-t4.*t24.*t27.*t42.*t43.*u5.*1.6e5).*(1.0./2.0e1)-t4.*t24.*t41+t4.*t24.*t56.*u5-t2.*t4.*t24.*t29.*u5.*2.0e1+t2.*t29.*t58.*t59.*u5.*2.0e1).*(1.0./1.0e1)-t8.*t16.*t32.*t35.*t48;-t8.*t16.*t32.*t35.*u9;-t8.*t16.*t32.*t35.*u10;-t8.*t16.*t32.*t35.*u11;t18.*t23.*(param8.*t6.*(t4.*t11.*t24.*t30.*t31.*u3.*1.6e4+t4.*t24.*t30.*t31.*t49.*u2.*1.6e4-t4.*t24.*t25.*t42.*t43.*u9.*1.6e5-t4.*t24.*t27.*t42.*t43.*u9.*1.6e5).*(1.0./2.0e1)-t4.*t24.*t52+t4.*t24.*t56.*u9-t2.*t4.*t24.*t29.*u9.*2.0e1+t2.*t29.*t58.*t59.*u9.*2.0e1).*(1.0./1.0e1)-t8.*t16.*t32.*t35.*t63;t8.*t16.*t21.*t22.*t32.*u5.^2.*-4.0e2;t8.*t16.*t21.*t22.*t32.*u5.*u6.*-4.0e2;t8.*t16.*t21.*t22.*t32.*u5.*u7.*-4.0e2;t18.*t23.*(param8.*t6.*(t9.*t21.*t22.*4.0e2-t30.*t31.*u2.*u5.*8.0e3).*(1.0./2.0e1)+t4.*t21.*t22.*t24.*u2.*u5.*4.0e2).*(-1.0./1.0e1)-t8.*t16.*t21.*t22.*t32.*t48.*u5.*4.0e2;t8.*t16.*t21.*t22.*t32.*u5.*u9.*-4.0e2;t8.*t16.*t21.*t22.*t32.*u5.*u10.*-4.0e2;t8.*t16.*t21.*t22.*t32.*u5.*u11.*-4.0e2;t18.*t23.*(param8.*t6.*(t21.*t22.*t49.*4.0e2-t30.*t31.*u2.*u9.*8.0e3).*(1.0./2.0e1)+t4.*t21.*t22.*t24.*u2.*u9.*4.0e2).*(-1.0./1.0e1)-t8.*t16.*t21.*t22.*t32.*t63.*u5.*4.0e2;t8.*t16.*t21.*t22.*t32.*u5.*u9.*-4.0e2;t8.*t16.*t21.*t22.*t32.*u6.*u9.*-4.0e2;t8.*t16.*t21.*t22.*t32.*u7.*u9.*-4.0e2;t18.*t23.*(param8.*t6.*(t21.*t22.*t39.*4.0e2-t30.*t31.*u3.*u5.*8.0e3).*(1.0./2.0e1)+t4.*t21.*t22.*t24.*u3.*u5.*4.0e2).*(-1.0./1.0e1)-t8.*t16.*t21.*t22.*t32.*t48.*u9.*4.0e2;t8.*t16.*t21.*t22.*t32.*u9.^2.*-4.0e2;t8.*t16.*t21.*t22.*t32.*u9.*u10.*-4.0e2;t8.*t16.*t21.*t22.*t32.*u9.*u11.*-4.0e2;t18.*t23.*(param8.*t6.*(t11.*t21.*t22.*4.0e2-t30.*t31.*u3.*u9.*8.0e3).*(1.0./2.0e1)+t4.*t21.*t22.*t24.*u3.*u9.*4.0e2).*(-1.0./1.0e1)-t8.*t16.*t21.*t22.*t32.*t63.*u9.*4.0e2;zero;zero;zero;zero;zero;zero;zero;zero;t66-t8.*t16.*t21.*t22.*t32.*u2.*u5.*4.0e2;t8.*t16.*t21.*t22.*t32.*u2.*u6.*-4.0e2;t8.*t16.*t21.*t22.*t32.*u2.*u7.*-4.0e2;t71-t8.*t16.*t21.*t22.*t32.*t48.*u2.*4.0e2;t8.*t16.*t21.*t22.*t32.*u2.*u9.*-4.0e2;t8.*t16.*t21.*t22.*t32.*u2.*u10.*-4.0e2;t8.*t16.*t21.*t22.*t32.*u2.*u11.*-4.0e2;t8.*t16.*t21.*t22.*t32.*t63.*u2.*-4.0e2;t72;t66+t73;t74;t75-t2.*t7.*t18.*t23.*u2.*2.0;t76;t77;t78;t79;zero;zero;t66;t2.*t7.*t18.*t23.*u3.*-2.0;zero;zero;zero;zero;zero;zero;zero;t80;zero;zero;zero;zero;t8.*t16.*t21.*t22.*t32.*u3.*u5.*-4.0e2;t8.*t16.*t21.*t22.*t32.*u3.*u6.*-4.0e2;t8.*t16.*t21.*t22.*t32.*u3.*u7.*-4.0e2;t8.*t16.*t21.*t22.*t32.*t48.*u3.*-4.0e2;t66-t8.*t16.*t21.*t22.*t32.*u3.*u9.*4.0e2;t8.*t16.*t21.*t22.*t32.*u3.*u10.*-4.0e2;t8.*t16.*t21.*t22.*t32.*u3.*u11.*-4.0e2;t71-t8.*t16.*t21.*t22.*t32.*t63.*u3.*4.0e2;zero;zero;zero;zero;zero;t66;zero;t2.*t7.*t18.*t23.*u2.*-2.0;t72;t73;t74;t75;t76;t77;t66+t78;t79-t2.*t7.*t18.*t23.*u3.*2.0;zero;zero;zero;zero;zero;zero;zero;t80];
end
f = reshape(f,ng,nch,nd);
f_udg = reshape(f_udg,ng,nch,nd,nc);
