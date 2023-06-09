function [f,f_udg] = flux_polarns2d(pg,udg,param,time)
%FLUX_POLARNS2D
%    [F,F_UDG] = FLUX_POLARNS2D(PG,UDG,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    12-Dec-2019 18:16:57
[ng,nc] = size(udg);
nch = 4;
nd = 2;
one = ones(ng,1);
param1 = param{1};
param3 = param{3};
param4 = param{4};
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
x1 = pg(:,1);
zero = zeros(ng,1);
t2 = 1.0./u1;
t3 = u2.^2;
t4 = 1.0./u1.^2;
t5 = 1.0./param3;
t6 = 1.0./x1;
t7 = t3.*t4.*(1.0./2.0);
t8 = u3.^2;
t9 = t4.*t8.*(1.0./2.0);
t10 = t7+t9;
t29 = t10.*u1;
t11 = -t29+u4;
t12 = param1-1.0;
t28 = t2.*u3.*u5;
t13 = -t28+u7;
t14 = t2.*t13;
t30 = t2.*u2.*u9;
t15 = -t30+u10;
t16 = t2.*t15;
t17 = t2.*u3;
t18 = t16+t17;
t19 = t6.*t18;
t20 = t14+t19;
t27 = t2.*u2.*u5;
t21 = -t27+u6;
t22 = t2.*t21.*2.0;
t34 = t2.*u3.*u9;
t23 = -t34+u11;
t24 = t2.*t23;
t50 = t2.*u2;
t25 = t24-t50;
t51 = t6.*t25;
t26 = t22-t51;
t31 = t5.*t20;
t32 = t2.*u2.*u3;
t33 = t31+t32;
t35 = t11.*t12;
t36 = t2.*u4;
t37 = t2.*t11.*t12;
t38 = t36+t37;
t39 = t2.*t21;
t40 = t2.*t23.*2.0;
t81 = t2.*u2.*2.0;
t41 = t40-t81;
t82 = t6.*t41;
t42 = t39-t82;
t43 = 1.0./param4;
t44 = 1.0./t12;
t45 = 1.0./u1.^3;
t46 = t3.*t45;
t47 = t8.*t45;
t48 = t46+t47;
t73 = t48.*u1;
t49 = t7+t9-t73;
t52 = t4.*t13;
t53 = t4.*t15;
t54 = t4.*u3;
t75 = t45.*u2.*u9;
t55 = t53+t54-t75;
t56 = t6.*t55;
t76 = t45.*u3.*u5;
t57 = t52+t56-t76;
t58 = t4.*t21.*2.0;
t59 = t4.*u2;
t60 = t45.*u3.*u9;
t109 = t4.*t23;
t61 = t59+t60-t109;
t62 = t6.*t61;
t63 = t58+t62-t45.*u2.*u5.*2.0;
t64 = t10.*u5;
t65 = t4.*t21.*u2;
t66 = t4.*t13.*u3;
t67 = t65+t66;
t68 = t67.*u1;
t69 = t64+t68-u8;
t70 = t12.*t69.*u1;
t71 = t11.*t12.*u5;
t72 = t70+t71;
t74 = 1.0./u1.^4;
t77 = t4.*u4;
t78 = t4.*t11.*t12;
t79 = t2.*t12.*t49;
t80 = t77+t78+t79;
t83 = t4.*t21;
t84 = t4.*u2.*2.0;
t85 = t45.*u3.*u9.*2.0;
t86 = t84+t85-t4.*t23.*2.0;
t87 = t6.*t86;
t100 = t45.*u2.*u5;
t88 = t83+t87-t100;
t89 = t10.*u9;
t90 = t4.*t15.*u2;
t91 = t4.*t23.*u3;
t92 = t90+t91;
t93 = t92.*u1;
t94 = t89+t93-u12;
t95 = t12.*t94.*u1;
t96 = t11.*t12.*u9;
t97 = t95+t96;
f = [u2;t35+t2.*t3+t5.*t26.*(2.0./3.0);t33;t38.*u2+t2.*t5.*t20.*u3+t2.*t5.*t26.*u2.*(2.0./3.0)-param1.*t4.*t5.*t43.*t44.*t72;t6.*u3;t6.*t33;t6.*(t35+t2.*t8-t5.*t42.*(2.0./3.0));t6.*(t38.*u3+t2.*t5.*t20.*u2-t2.*t5.*t42.*u3.*(2.0./3.0)-param1.*t4.*t5.*t6.*t43.*t44.*t97)];
if nargout > 1
    t98 = t2.*t6;
    t99 = t98-t4.*u5.*2.0;
    t101 = t17-t4.*t5.*t6.*u9;
    t102 = t2.*t6.*2.0;
    t104 = t4.*u5;
    t103 = t102-t104;
    t105 = t98-t104;
    t106 = t4.*t12.*u2.*u3;
    t107 = t5.*t105;
    t108 = t50+t107;
    t110 = t2.*t12;
    t111 = t2+t110;
    t112 = t4.*t5.*t6.*u2;
    t113 = 1.0./x1.^2;
    t114 = t12.*t49.*u1;
    t115 = t35+t114;
    t116 = t2.*t5.*t6;
    t117 = t4.*t5.*t6.*u2.*(2.0./3.0);
    f_udg = [zero;-t3.*t4-t12.*t49-t5.*t63.*(2.0./3.0);-t5.*t57-t4.*u2.*u3;-t80.*u2-t4.*t5.*t20.*u3-t4.*t5.*t26.*u2.*(2.0./3.0)-t2.*t5.*t57.*u3-t2.*t5.*t63.*u2.*(2.0./3.0)-param1.*t4.*t5.*t43.*t44.*(t12.*t69-t12.*t49.*u5+t12.*u1.*(t65+t66-u1.*(t13.*t45.*u3.*2.0+t21.*t45.*u2.*2.0-t3.*t74.*u5-t8.*t74.*u5)-t48.*u5))+param1.*t5.*t43.*t44.*t45.*t72.*2.0;zero;-t6.*(t5.*t57+t4.*u2.*u3);-t6.*(t4.*t8+t12.*t49-t5.*t88.*(2.0./3.0));-t6.*(t80.*u3+t4.*t5.*t20.*u2-t4.*t5.*t42.*u3.*(2.0./3.0)+t2.*t5.*t57.*u2-t2.*t5.*t88.*u3.*(2.0./3.0)+param1.*t4.*t5.*t6.*t43.*t44.*(t12.*t94-t12.*t49.*u9+t12.*u1.*(t90+t91-u1.*(t15.*t45.*u2.*2.0+t23.*t45.*u3.*2.0-t3.*t74.*u9-t8.*t74.*u9)-t48.*u9))-param1.*t5.*t6.*t43.*t44.*t45.*t97.*2.0);one;t81+t5.*t99.*(2.0./3.0)-t2.*t12.*u2;t101;t36+t37-t3.*t4.*t12+t2.*t5.*t26.*(2.0./3.0)+t2.*t5.*t99.*u2.*(2.0./3.0)-t5.*t6.*t45.*u3.*u9-param1.*t4.*t5.*t43.*t44.*(t12.*u1.*(u1.*(t83-t100)+t4.*u2.*u5)-t2.*t12.*u2.*u5);zero;t6.*t101;-t6.*(t5.*t103.*(2.0./3.0)+t2.*t12.*u2);-t6.*(t106-t2.*t5.*t20+t2.*t5.*t103.*u3.*(2.0./3.0)+t5.*t6.*t45.*u2.*u9+param1.*t4.*t5.*t6.*t43.*t44.*(t12.*u1.*(u1.*(t53-t75)+t4.*u2.*u9)-t2.*t12.*u2.*u9));zero;-t2.*t12.*u3+t4.*t5.*t6.*u9.*(2.0./3.0);t108;-t106+t2.*t5.*t20+t2.*t5.*t105.*u3+t5.*t6.*t45.*u2.*u9.*(2.0./3.0)-param1.*t4.*t5.*t43.*t44.*(t12.*u1.*(u1.*(t52-t76)+t4.*u3.*u5)-t2.*t12.*u3.*u5);one.*t6;t6.*t108;-t6.*(t2.*u3.*-2.0+t2.*t12.*u3+t4.*t5.*t6.*u9.*(4.0./3.0));t6.*(t36+t37-t4.*t8.*t12-t2.*t5.*t42.*(2.0./3.0)+t2.*t5.*t105.*u2-t5.*t6.*t45.*u3.*u9.*(4.0./3.0)+param1.*t4.*t5.*t6.*t43.*t44.*(t12.*u1.*(u1.*(t60-t109)-t4.*u3.*u9)+t2.*t12.*u3.*u9));zero;one.*t12;zero;t111.*u2-param1.*t4.*t5.*t43.*u5;zero;zero;one.*t6.*t12;t6.*(t111.*u3-param1.*t4.*t5.*t6.*t43.*u9);zero;t4.*t5.*u2.*(-4.0./3.0);-t4.*t5.*u3;t3.*t5.*t45.*(-4.0./3.0)-t5.*t8.*t45-param1.*t4.*t5.*t43.*t44.*t115;zero;-t4.*t5.*t6.*u3;t117;t5.*t6.*t45.*u2.*u3.*(-1.0./3.0);zero;t2.*t5.*(4.0./3.0);zero;t4.*t5.*u2.*(4.0./3.0)-param1.*t4.*t5.*t43.*u2;zero;zero;t2.*t5.*t6.*(-2.0./3.0);t4.*t5.*t6.*u3.*(-2.0./3.0);zero;zero;t2.*t5;t4.*t5.*u3-param1.*t4.*t5.*t43.*u3;zero;t116;zero;t112;zero;zero;zero;param1.*t2.*t5.*t43;zero;zero;zero;zero;zero;t4.*t5.*t6.*u3.*(2.0./3.0);-t112;t5.*t6.*t45.*u2.*u3.*(-1.0./3.0);zero;-t4.*t5.*t113.*u2;t4.*t5.*t113.*u3.*(-4.0./3.0);-t6.*(t3.*t5.*t6.*t45+t5.*t6.*t8.*t45.*(4.0./3.0)+param1.*t4.*t5.*t6.*t43.*t44.*t115);zero;zero;t116;t4.*t5.*t6.*u3;zero;t2.*t5.*t113;zero;t6.*(t112-param1.*t4.*t5.*t6.*t43.*u2);zero;t2.*t5.*t6.*(-2.0./3.0);zero;-t117;zero;zero;t2.*t5.*t113.*(4.0./3.0);t6.*(t4.*t5.*t6.*u3.*(4.0./3.0)-param1.*t4.*t5.*t6.*t43.*u3);zero;zero;zero;zero;zero;zero;zero;param1.*t2.*t5.*t43.*t113];
end
f = reshape(f,ng,nch,nd);
f_udg = reshape(f_udg,ng,nch,nd,nc);
