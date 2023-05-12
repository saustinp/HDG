function [f,f_udg] = avflux2d4(pg,udg,param,time)
%AVFLUX2D4
%    [F,F_UDG] = AVFLUX2D4(PG,UDG,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    27-Mar-2013 09:10:44
[ng,nc] = size(udg);
nch = 4;
nd = 2;
param1 = param{1};
param8 = param{8};
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
t2 = 1.0./param8;
t3 = t2.*u1.*2.0e1;
t4 = exp(t3);
t5 = t4+2.0e1;
t6 = log(t5);
t7 = 1.0./t6;
t8 = 1.0./param8.^2;
t9 = 1.0./t6.^2;
t10 = u2.^2;
t11 = u3.^2;
t12 = param1-1.0;
t37 = t2.*t7.*u2.*u5.*2.0e1;
t13 = -t37+u6;
t14 = t2.*t7.*t13.*2.0e1;
t38 = t2.*t7.*u3.*u9.*2.0e1;
t15 = -t38+u11;
t16 = t2.*t7.*t15.*2.0e1;
t17 = t14+t16;
t18 = t8.*t9.*t10.*4.0e2;
t19 = t8.*t9.*t11.*4.0e2;
t20 = t18+t19;
t21 = t12.*t20;
t22 = t8.*t9.*t10.*2.0e2;
t23 = t8.*t9.*t11.*2.0e2;
t24 = t22+t23;
t39 = param8.*t6.*t24.*(1.0./2.0e1);
t25 = -t39+u4;
t26 = param1.*t2.*t7.*t12.*t25.*4.0e1;
t27 = t21+t26;
t28 = param1+1.0;
t29 = 1.0./t28;
t30 = t27.*t29;
t31 = 1.0./sqrt(t30);
t32 = t17.*t31.*1.0e1;
t33 = t32-5.0;
t34 = exp(t33);
t35 = t34+1.0;
t36 = log(t35);
t40 = 1.0./t5;
t41 = 1.0./param8.^3;
t42 = 1.0./t6.^3;
t43 = 1.0./t35;
t44 = t4.*t8.*t9.*t13.*t40.*4.0e2;
t45 = t4.*t8.*t9.*t15.*t40.*4.0e2;
t62 = t4.*t40.*t41.*t42.*u2.*u5.*8.0e3;
t63 = t4.*t40.*t41.*t42.*u3.*u9.*8.0e3;
t46 = t44+t45-t62-t63;
t47 = t31.*t46.*1.0e1;
t48 = 1.0./t30.^(3.0./2.0);
t49 = t4.*t10.*t40.*t41.*t42.*1.6e4;
t50 = t4.*t11.*t40.*t41.*t42.*1.6e4;
t51 = t49+t50;
t52 = t12.*t51;
t53 = t4.*t24.*t40;
t54 = t4.*t10.*t40.*t41.*t42.*8.0e3;
t55 = t4.*t11.*t40.*t41.*t42.*8.0e3;
t56 = t54+t55;
t64 = param8.*t6.*t56.*(1.0./2.0e1);
t57 = t53-t64;
t58 = param1.*t2.*t7.*t12.*t57.*4.0e1;
t59 = param1.*t4.*t8.*t9.*t12.*t25.*t40.*8.0e2;
t60 = t52+t58+t59;
t65 = t17.*t29.*t48.*t60.*5.0;
t61 = t47-t65;
t66 = t8.*t9.*t13.*u2.*4.0e2;
t72 = t2.*t7.*u3.*u5.*2.0e1;
t67 = -t72+u7;
t68 = t8.*t9.*t67.*u3.*4.0e2;
t69 = t66+t68;
t70 = 1.0./param8.^4;
t71 = 1.0./t6.^4;
t73 = param8.*t6.*t69.*(1.0./2.0e1);
t74 = t4.*t24.*t40.*u5;
t75 = t73+t74-u8;
t95 = t12.*t75;
t76 = -t95+u8;
t81 = t2.*t7.*u2.*u9.*2.0e1;
t77 = -t81+u10;
t78 = t8.*t9.*t77.*u2.*4.0e2;
t79 = t8.*t9.*t15.*u3.*4.0e2;
t80 = t78+t79;
t82 = t2.*u1.*4.0e1;
t83 = exp(t82);
t84 = 1.0./t5.^2;
t85 = param8.*t6.*t80.*(1.0./2.0e1);
t86 = t4.*t24.*t40.*u9;
t87 = t85+t86-u12;
t96 = t12.*t87;
t88 = -t96+u12;
f = [t36.*u5.*(1.0./1.0e1);t36.*u6.*(1.0./1.0e1);t36.*u7.*(1.0./1.0e1);t36.*t76.*(1.0./1.0e1);t36.*u9.*(1.0./1.0e1);t36.*u10.*(1.0./1.0e1);t36.*u11.*(1.0./1.0e1);t36.*t88.*(1.0./1.0e1)];
if nargout > 1
    t89 = t8.*t9.*t31.*u5.*4.0e3;
    t90 = t8.*t9.*t12.*u2.*8.0e2;
    t94 = param1.*t8.*t9.*t12.*u2.*8.0e2;
    t91 = t90-t94;
    t92 = t17.*t29.*t48.*t91.*5.0;
    t93 = t89+t92;
    t97 = t8.*t9.*t31.*u9.*4.0e3;
    t98 = t8.*t9.*t12.*u3.*8.0e2;
    t102 = param1.*t8.*t9.*t12.*u3.*8.0e2;
    t99 = t98-t102;
    t100 = t17.*t29.*t48.*t99.*5.0;
    t101 = t97+t100;
    t103 = t36.*(1.0./1.0e1);
    t104 = t10.*t41.*t42.*8.0e3;
    t105 = t11.*t41.*t42.*8.0e3;
    t106 = t104+t105;
    t107 = t53-param8.*t6.*t106.*(1.0./2.0e1);
    t108 = t2.*t7.*t31.*t34.*t43.*u5.*2.0e1;
    t109 = t2.*t7.*t31.*t34.*t43.*u6.*2.0e1;
    t110 = t2.*t7.*t31.*t34.*t43.*u7.*2.0e1;
    t111 = t2.*t7.*t31.*t34.*t43.*t76.*2.0e1;
    t112 = t2.*t7.*t31.*t34.*t43.*u9.*2.0e1;
    t113 = t2.*t7.*t31.*t34.*t43.*u10.*2.0e1;
    t114 = t2.*t7.*t31.*t34.*t43.*u11.*2.0e1;
    t115 = t2.*t7.*t31.*t34.*t43.*t88.*2.0e1;
    t116 = param1.*t36.*(1.0./1.0e1);
    f_udg = [t34.*t43.*t61.*u5.*(-1.0./1.0e1);t34.*t43.*t61.*u6.*(-1.0./1.0e1);t34.*t43.*t61.*u7.*(-1.0./1.0e1);t12.*t36.*(param8.*t6.*(t4.*t13.*t40.*t41.*t42.*u2.*1.6e4+t4.*t40.*t41.*t42.*t67.*u3.*1.6e4-t4.*t10.*t40.*t70.*t71.*u5.*1.6e5-t4.*t11.*t40.*t70.*t71.*u5.*1.6e5).*(1.0./2.0e1)-t4.*t40.*t69+t4.*t40.*t56.*u5-t2.*t4.*t24.*t40.*u5.*2.0e1+t2.*t24.*t83.*t84.*u5.*2.0e1).*(1.0./1.0e1)-t34.*t43.*t61.*t76.*(1.0./1.0e1);t34.*t43.*t61.*u9.*(-1.0./1.0e1);t34.*t43.*t61.*u10.*(-1.0./1.0e1);t34.*t43.*t61.*u11.*(-1.0./1.0e1);t12.*t36.*(param8.*t6.*(t4.*t15.*t40.*t41.*t42.*u3.*1.6e4-t4.*t10.*t40.*t70.*t71.*u9.*1.6e5-t4.*t11.*t40.*t70.*t71.*u9.*1.6e5+t4.*t40.*t41.*t42.*t77.*u2.*1.6e4).*(1.0./2.0e1)-t4.*t40.*t80+t4.*t40.*t56.*u9-t2.*t4.*t24.*t40.*u9.*2.0e1+t2.*t24.*t83.*t84.*u9.*2.0e1).*(1.0./1.0e1)-t34.*t43.*t61.*t88.*(1.0./1.0e1);t34.*t43.*t93.*u5.*(-1.0./1.0e1);t34.*t43.*t93.*u6.*(-1.0./1.0e1);t34.*t43.*t93.*u7.*(-1.0./1.0e1);t12.*t36.*(param8.*t6.*(t8.*t9.*t13.*4.0e2-t41.*t42.*u2.*u5.*8.0e3).*(1.0./2.0e1)+t4.*t8.*t9.*t40.*u2.*u5.*4.0e2).*(-1.0./1.0e1)-t34.*t43.*t76.*t93.*(1.0./1.0e1);t34.*t43.*t93.*u9.*(-1.0./1.0e1);t34.*t43.*t93.*u10.*(-1.0./1.0e1);t34.*t43.*t93.*u11.*(-1.0./1.0e1);t12.*t36.*(param8.*t6.*(t8.*t9.*t77.*4.0e2-t41.*t42.*u2.*u9.*8.0e3).*(1.0./2.0e1)+t4.*t8.*t9.*t40.*u2.*u9.*4.0e2).*(-1.0./1.0e1)-t34.*t43.*t88.*t93.*(1.0./1.0e1);t34.*t43.*t101.*u5.*(-1.0./1.0e1);t34.*t43.*t101.*u6.*(-1.0./1.0e1);t34.*t43.*t101.*u7.*(-1.0./1.0e1);t12.*t36.*(param8.*t6.*(t8.*t9.*t67.*4.0e2-t41.*t42.*u3.*u5.*8.0e3).*(1.0./2.0e1)+t4.*t8.*t9.*t40.*u3.*u5.*4.0e2).*(-1.0./1.0e1)-t34.*t43.*t76.*t101.*(1.0./1.0e1);t34.*t43.*t101.*u9.*(-1.0./1.0e1);t34.*t43.*t101.*u10.*(-1.0./1.0e1);t34.*t43.*t101.*u11.*(-1.0./1.0e1);t12.*t36.*(param8.*t6.*(t8.*t9.*t15.*4.0e2-t41.*t42.*u3.*u9.*8.0e3).*(1.0./2.0e1)+t4.*t8.*t9.*t40.*u3.*u9.*4.0e2).*(-1.0./1.0e1)-t34.*t43.*t88.*t101.*(1.0./1.0e1);param1.*t2.*t7.*t12.*t17.*t29.*t34.*t43.*t48.*u5.*-2.0e1;param1.*t2.*t7.*t12.*t17.*t29.*t34.*t43.*t48.*u6.*-2.0e1;param1.*t2.*t7.*t12.*t17.*t29.*t34.*t43.*t48.*u7.*-2.0e1;param1.*t2.*t7.*t12.*t17.*t29.*t34.*t43.*t48.*t76.*-2.0e1;param1.*t2.*t7.*t12.*t17.*t29.*t34.*t43.*t48.*u9.*-2.0e1;param1.*t2.*t7.*t12.*t17.*t29.*t34.*t43.*t48.*u10.*-2.0e1;param1.*t2.*t7.*t12.*t17.*t29.*t34.*t43.*t48.*u11.*-2.0e1;param1.*t2.*t7.*t12.*t17.*t29.*t34.*t43.*t48.*t88.*-2.0e1;t103-t8.*t9.*t31.*t34.*t43.*u2.*u5.*4.0e2;t8.*t9.*t31.*t34.*t43.*u2.*u6.*-4.0e2;t8.*t9.*t31.*t34.*t43.*u2.*u7.*-4.0e2;t12.*t36.*t107.*(-1.0./1.0e1)-t8.*t9.*t31.*t34.*t43.*t76.*u2.*4.0e2;t8.*t9.*t31.*t34.*t43.*u2.*u9.*-4.0e2;t8.*t9.*t31.*t34.*t43.*u2.*u10.*-4.0e2;t8.*t9.*t31.*t34.*t43.*u2.*u11.*-4.0e2;t8.*t9.*t31.*t34.*t43.*t88.*u2.*-4.0e2;t108;t103+t109;t110;t111-t2.*t7.*t12.*t36.*u2.*2.0;t112;t113;t114;t115;zero;zero;t103;t2.*t7.*t12.*t36.*u3.*-2.0;zero;zero;zero;zero;zero;zero;zero;t116;zero;zero;zero;zero;t8.*t9.*t31.*t34.*t43.*u3.*u5.*-4.0e2;t8.*t9.*t31.*t34.*t43.*u3.*u6.*-4.0e2;t8.*t9.*t31.*t34.*t43.*u3.*u7.*-4.0e2;t8.*t9.*t31.*t34.*t43.*t76.*u3.*-4.0e2;t103-t8.*t9.*t31.*t34.*t43.*u3.*u9.*4.0e2;t8.*t9.*t31.*t34.*t43.*u3.*u10.*-4.0e2;t8.*t9.*t31.*t34.*t43.*u3.*u11.*-4.0e2;t12.*t36.*t107.*(-1.0./1.0e1)-t8.*t9.*t31.*t34.*t43.*t88.*u3.*4.0e2;zero;zero;zero;zero;zero;t103;zero;t2.*t7.*t12.*t36.*u2.*-2.0;t108;t109;t110;t111;t112;t113;t103+t114;t115-t2.*t7.*t12.*t36.*u3.*2.0;zero;zero;zero;zero;zero;zero;zero;t116];
end
f = reshape(f,ng,nch,nd);
f_udg = reshape(f_udg,ng,nch,nd,nc);
