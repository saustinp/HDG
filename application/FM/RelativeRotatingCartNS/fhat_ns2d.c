function [fh,fh_udg,fh_uh] = fhat_ns2d.c(nl,pg,udg,uh,param,time)
%FHAT_NS2D.C
%    [FH,FH_UDG,FH_UH] = FHAT_NS2D.C(NL,PG,UDG,UH,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    13-Dec-2019 21:40:22
[ng,nc] = size(udg);
nch = 4;
nd = 2;
one = ones(ng,1);
param1 = param{1};
param3 = param{3};
param4 = param{4};
u5 = udg(:,5);
u6 = udg(:,6);
u7 = udg(:,7);
u8 = udg(:,8);
u9 = udg(:,9);
u10 = udg(:,10);
u11 = udg(:,11);
u12 = udg(:,12);
uh1 = uh(:,1);
uh2 = uh(:,2);
uh3 = uh(:,3);
uh4 = uh(:,4);
zero = zeros(ng,1);
t2 = one.^2;
t3 = t2.^2;
t4 = uh2.^2;
t5 = 1.0./uh1.^2;
t6 = 1.0./uh1;
t7 = 1.0./param3;
t8 = t4.*t5.*(1.0./2.0);
t9 = uh3.^2;
t10 = t5.*t9.*(1.0./2.0);
t11 = t8+t10;
t25 = t11.*uh1;
t12 = -t25+uh4;
t13 = param1-1.0;
t24 = t6.*u5.*uh3;
t14 = -t24+u7;
t15 = t6.*t14;
t26 = t6.*u9.*uh2;
t16 = -t26+u10;
t17 = t6.*t16;
t18 = t15+t17;
t23 = t6.*u5.*uh2;
t19 = -t23+u6;
t20 = t6.*t19.*2.0;
t31 = t6.*u9.*uh3;
t21 = -t31+u11;
t53 = t6.*t21;
t22 = t20-t53;
t27 = t7.*t18;
t28 = t6.*uh2.*uh3;
t29 = t27+t28;
t30 = t12.*t13;
t32 = t6.*uh4;
t33 = t6.*t12.*t13;
t34 = t32+t33;
t35 = t6.*t19;
t77 = t6.*t21.*2.0;
t36 = t35-t77;
t37 = 1.0./param4;
t38 = 1.0./t13;
t39 = 1.0./uh1.^3;
t40 = one.*t6.*t7;
t41 = one.*t5.*t7.*uh2;
t42 = t4.*t39;
t43 = t9.*t39;
t44 = t42+t43;
t52 = t44.*uh1;
t45 = t8+t10-t52;
t46 = t13.*t45.*uh1;
t47 = t30+t46;
t48 = param1.*t5.*t7.*t37.*t38.*t47;
t49 = one.*t5.*t7.*uh2.*(2.0./3.0);
t50 = one.*t6.*t7.*(4.0./3.0);
t51 = one.*param1.*t6.*t7.*t37;
t54 = t5.*t14;
t55 = t5.*t16;
t70 = t39.*u5.*uh3;
t71 = t39.*u9.*uh2;
t56 = t54+t55-t70-t71;
t57 = t5.*t19.*2.0;
t58 = t39.*u9.*uh3;
t96 = t5.*t21;
t59 = t57+t58-t96-t39.*u5.*uh2.*2.0;
t60 = t11.*u5;
t61 = t5.*t19.*uh2;
t62 = t5.*t14.*uh3;
t63 = t61+t62;
t64 = t63.*uh1;
t65 = t60+t64-u8;
t66 = t13.*t65.*uh1;
t67 = t12.*t13.*u5;
t68 = t66+t67;
t69 = 1.0./uh1.^4;
t72 = -t7.*t56-t5.*uh2.*uh3;
t73 = t5.*uh4;
t74 = t5.*t12.*t13;
t75 = t6.*t13.*t45;
t76 = t73+t74+t75;
t78 = t5.*t19;
t79 = t39.*u9.*uh3.*2.0;
t90 = t39.*u5.*uh2;
t80 = t78+t79-t90-t5.*t21.*2.0;
t81 = t11.*u9;
t82 = t5.*t16.*uh2;
t83 = t5.*t21.*uh3;
t84 = t82+t83;
t85 = t84.*uh1;
t86 = t81+t85-u12;
t87 = t13.*t86.*uh1;
t88 = t12.*t13.*u9;
t89 = t87+t88;
fh = [t2.*t3.*uh2;t30+t4.*t6+t7.*t22.*(2.0./3.0);t29;t34.*uh2+t6.*t7.*t18.*uh3+t6.*t7.*t22.*uh2.*(2.0./3.0)-param1.*t5.*t7.*t37.*t38.*t68;t2.*t3.*uh3;t29;t30+t6.*t9-t7.*t36.*(2.0./3.0);t34.*uh3+t6.*t7.*t18.*uh2-t6.*t7.*t36.*uh3.*(2.0./3.0)-param1.*t5.*t7.*t37.*t38.*t89];
if nargout > 1
    fh_udg = [zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;one.*t5.*t7.*uh2.*(-4.0./3.0);-one.*t5.*t7.*uh3;-one.*(t48+t4.*t7.*t39.*(4.0./3.0)+t7.*t9.*t39);zero;-one.*t5.*t7.*uh3;t49;one.*t7.*t39.*uh2.*uh3.*(-1.0./3.0);zero;t50;zero;one.*(t5.*t7.*uh2.*(4.0./3.0)-param1.*t5.*t7.*t37.*uh2);zero;zero;one.*t6.*t7.*(-2.0./3.0);one.*t5.*t7.*uh3.*(-2.0./3.0);zero;zero;t40;one.*(t5.*t7.*uh3-param1.*t5.*t7.*t37.*uh3);zero;t40;zero;t41;zero;zero;zero;t51;zero;zero;zero;zero;zero;one.*t5.*t7.*uh3.*(2.0./3.0);-t41;one.*t7.*t39.*uh2.*uh3.*(-1.0./3.0);zero;-t41;one.*t5.*t7.*uh3.*(-4.0./3.0);-one.*(t48+t4.*t7.*t39+t7.*t9.*t39.*(4.0./3.0));zero;zero;t40;one.*t5.*t7.*uh3;zero;t40;zero;one.*(t5.*t7.*uh2-param1.*t5.*t7.*t37.*uh2);zero;one.*t6.*t7.*(-2.0./3.0);zero;-t49;zero;zero;t50;one.*(t5.*t7.*uh3.*(4.0./3.0)-param1.*t5.*t7.*t37.*uh3);zero;zero;zero;zero;zero;zero;zero;t51];
end
if nargout > 2
    t91 = t6.*uh3;
    t92 = t91-t5.*t7.*u9;
    t93 = t6.*t7.*t18;
    t94 = t6.*uh2;
    t95 = t94-t5.*t7.*u5;
    t97 = one.*t13;
    t98 = t6.*t13;
    t99 = t6+t98;
    fh_uh = [zero;-t4.*t5-t13.*t45-t7.*t59.*(2.0./3.0);t72;-t76.*uh2-t5.*t7.*t18.*uh3-t5.*t7.*t22.*uh2.*(2.0./3.0)-t6.*t7.*t56.*uh3-t6.*t7.*t59.*uh2.*(2.0./3.0)-param1.*t5.*t7.*t37.*t38.*(t13.*t65-t13.*t45.*u5+t13.*uh1.*(t61+t62+uh1.*(t4.*t69.*u5+t9.*t69.*u5-t14.*t39.*uh3.*2.0-t19.*t39.*uh2.*2.0)-t44.*u5))+param1.*t7.*t37.*t38.*t39.*t68.*2.0;zero;t72;-t5.*t9-t13.*t45+t7.*t80.*(2.0./3.0);-t76.*uh3-t5.*t7.*t18.*uh2+t5.*t7.*t36.*uh3.*(2.0./3.0)-t6.*t7.*t56.*uh2+t6.*t7.*t80.*uh3.*(2.0./3.0)-param1.*t5.*t7.*t37.*t38.*(t13.*t86-t13.*t45.*u9+t13.*uh1.*(t82+t83+uh1.*(t4.*t69.*u9+t9.*t69.*u9-t16.*t39.*uh2.*2.0-t21.*t39.*uh3.*2.0)-t44.*u9))+param1.*t7.*t37.*t38.*t39.*t89.*2.0;t2.*t3;t6.*uh2.*2.0-t5.*t7.*u5.*(4.0./3.0)-t6.*t13.*uh2;t92;t32+t33-t4.*t5.*t13+t6.*t7.*t22.*(2.0./3.0)-t7.*t39.*u5.*uh2.*(4.0./3.0)-t7.*t39.*u9.*uh3-param1.*t5.*t7.*t37.*t38.*(t13.*uh1.*(uh1.*(t78-t90)+t5.*u5.*uh2)-t6.*t13.*u5.*uh2);zero;t92;t5.*t7.*u5.*(2.0./3.0)-t6.*t13.*uh2;t93+t7.*t39.*u5.*uh3.*(2.0./3.0)-t7.*t39.*u9.*uh2-t5.*t13.*uh2.*uh3-param1.*t5.*t7.*t37.*t38.*(t13.*uh1.*(uh1.*(t55-t71)+t5.*u9.*uh2)-t6.*t13.*u9.*uh2);zero;t5.*t7.*u9.*(2.0./3.0)-t6.*t13.*uh3;t95;t93-t7.*t39.*u5.*uh3+t7.*t39.*u9.*uh2.*(2.0./3.0)-t5.*t13.*uh2.*uh3-param1.*t5.*t7.*t37.*t38.*(t13.*uh1.*(uh1.*(t54-t70)+t5.*u5.*uh3)-t6.*t13.*u5.*uh3);t2.*t3;t95;t6.*uh3.*2.0-t5.*t7.*u9.*(4.0./3.0)-t6.*t13.*uh3;t32+t33-t5.*t9.*t13-t6.*t7.*t36.*(2.0./3.0)-t7.*t39.*u5.*uh2-t7.*t39.*u9.*uh3.*(4.0./3.0)+param1.*t5.*t7.*t37.*t38.*(t13.*uh1.*(uh1.*(t58-t96)-t5.*u9.*uh3)+t6.*t13.*u9.*uh3);zero;t97;zero;t99.*uh2-param1.*t5.*t7.*t37.*u5;zero;zero;t97;t99.*uh3-param1.*t5.*t7.*t37.*u9];
end
fh = reshape(fh,ng,nch);
fh_udg = reshape(fh_udg,ng,nch,nc);
fh_uh = reshape(fh_uh,ng,nch,nch);
