function [s,s_udg] = source_ns3d(pg,udg,param,time)
%SOURCE_NS3D
%    [S,S_UDG] = SOURCE_NS3D(PG,UDG,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    15-Dec-2019 23:46:05
[ng,nc] = size(udg);
nch = 5;
nd = 3;
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
u13 = udg(:,13);
u16 = udg(:,16);
u17 = udg(:,17);
u19 = udg(:,19);
x1 = pg(:,1);
zero = zeros(ng,1);
t2 = 1.0./x1;
t3 = 1.0./u1;
t4 = 1.0./param3;
t24 = t3.*u2.*u6;
t5 = -t24+u7;
t25 = t3.*u4.*u16;
t6 = -t25+u19;
t27 = t3.*u3.*u11;
t7 = -t27+u13;
t8 = u2.^2;
t9 = 1.0./u1.^2;
t10 = u3.^2;
t34 = t3.*u4.*u6;
t11 = -t34+u9;
t12 = t3.*t11;
t44 = t3.*u2.*u16;
t13 = -t44+u17;
t14 = t3.*t13;
t15 = t12+t14;
t33 = t3.*u3.*u6;
t16 = -t33+u8;
t17 = t3.*t16;
t43 = t3.*u2.*u11;
t18 = -t43+u12;
t19 = t3.*t18;
t20 = t3.*u3;
t21 = t19+t20;
t22 = t2.*t21;
t23 = t17+t22;
t26 = t3.*t6;
t28 = t3.*t7;
t50 = t3.*u2;
t29 = t28-t50;
t30 = t2.*t29;
t49 = t3.*t5.*2.0;
t31 = t26+t30-t49;
t32 = param1-1.0;
t35 = t8.*t9.*(1.0./2.0);
t36 = t9.*t10.*(1.0./2.0);
t37 = u4.^2;
t38 = t9.*t37.*(1.0./2.0);
t39 = t35+t36+t38;
t45 = t39.*u1;
t40 = -t45+u5;
t41 = 1.0./u1.^3;
t42 = t41.*u4.*u16;
t46 = t9.*t11;
t47 = t9.*t13;
t91 = t41.*u4.*u6;
t48 = t46+t47-t91-t41.*u2.*u16;
t51 = t9.*t16;
t52 = t9.*t18;
t53 = t9.*u3;
t54 = t52+t53-t41.*u2.*u11;
t55 = t2.*t54;
t90 = t41.*u3.*u6;
t56 = t51+t55-t90;
t57 = t5.*t9.*2.0;
t58 = t6.*t9;
t59 = t9.*u2;
t60 = t41.*u3.*u11;
t61 = t59+t60-t7.*t9;
t62 = t2.*t61;
t63 = 1.0./param4;
t64 = 1.0./t32;
t65 = t5.*t9.*u2;
t66 = t9.*t16.*u3;
t67 = t9.*t11.*u4;
t68 = t65+t66+t67;
t69 = t68.*u1;
t70 = t39.*u6;
t71 = t69+t70-u10;
t72 = t8.*t41;
t73 = t10.*t41;
t74 = t37.*t41;
t75 = t72+t73+t74;
t92 = t75.*u1;
t76 = t35+t36+t38-t92;
t77 = 1.0./u1.^4;
t78 = t32.*t71.*u1;
t79 = t32.*t40.*u6;
t80 = t78+t79;
t81 = t3.*u5;
t82 = t3.*t32.*t40;
s = [-t2.*u2;-t2.*(t3.*t8-t3.*t10-t4.*t31.*(2.0./3.0)+t4.*(t26+t3.*t5-t2.*(t3.*t7.*2.0-t3.*u2.*2.0)).*(2.0./3.0));-t2.*(t4.*t23.*2.0+t3.*u2.*u3.*2.0);-t2.*(t4.*t15+t3.*u2.*u4);-t2.*(u2.*(t81+t82)+t3.*t4.*t15.*u4+t3.*t4.*t23.*u3-t3.*t4.*t31.*u2.*(2.0./3.0)-param1.*t4.*t9.*t63.*t64.*t80)];
if nargout > 1
    t83 = t2.*t3;
    t84 = t83-t9.*u6.*2.0;
    t85 = t5.*t9;
    t86 = t3.*u3.*2.0;
    t87 = t86-t2.*t4.*t9.*u11.*2.0;
    t88 = t3.*u2.*2.0;
    t89 = t83-t9.*u6;
    t93 = 1.0./x1.^2;
    t94 = t2.*t4.*t9.*u4;
    s_udg = [zero;t2.*(t4.*(-t42+t58+t85+t2.*(t7.*t9.*-2.0+t9.*u2.*2.0+t41.*u3.*u11.*2.0)-t41.*u2.*u6).*(2.0./3.0)+t8.*t9-t9.*t10+t4.*(t42+t57+t62-t6.*t9-t41.*u2.*u6.*2.0).*(2.0./3.0));t2.*(t4.*t56.*2.0+t9.*u2.*u3.*2.0);t2.*(t4.*t48+t9.*u2.*u4);t2.*(u2.*(t9.*u5+t9.*t32.*t40+t3.*t32.*t76)+t3.*t4.*u2.*(t42+t57-t58+t62-t41.*u2.*u6.*2.0).*(2.0./3.0)+t4.*t9.*t15.*u4+t4.*t9.*t23.*u3-t4.*t9.*t31.*u2.*(2.0./3.0)+t3.*t4.*t48.*u4+t3.*t4.*t56.*u3-param1.*t4.*t41.*t63.*t64.*t80.*2.0+param1.*t4.*t9.*t63.*t64.*(t32.*t71-t32.*t76.*u6+t32.*u1.*(t65+t66+t67-t75.*u6-u1.*(t5.*t41.*u2.*2.0+t11.*t41.*u4.*2.0+t16.*t41.*u3.*2.0-t8.*t77.*u6-t10.*t77.*u6-t37.*t77.*u6))));-one.*t2;-t2.*(t88+t4.*t84.*(2.0./3.0)+t4.*(t2.*t3.*2.0-t9.*u6).*(2.0./3.0));-t2.*t87;-t2.*(t3.*u4-t4.*t9.*u16);t2.*(-t81-t82+t3.*t4.*t31.*(2.0./3.0)+t8.*t9.*t32-t3.*t4.*t84.*u2.*(2.0./3.0)+t4.*t41.*u4.*u16+t2.*t4.*t41.*u3.*u11+param1.*t4.*t9.*t63.*t64.*(t32.*u1.*(u1.*(t85-t41.*u2.*u6)+t9.*u2.*u6)-t3.*t32.*u2.*u6));zero;t2.*t87;-t2.*(t88+t4.*t89.*2.0);zero;-t2.*(t3.*t4.*t23+t3.*t4.*t89.*u3-t9.*t32.*u2.*u3+t2.*t4.*t41.*u2.*u11.*(2.0./3.0)-param1.*t4.*t9.*t63.*t64.*(t32.*u1.*(u1.*(t51-t90)+t9.*u3.*u6)-t3.*t32.*u3.*u6));zero;zero;zero;-t2.*(t50-t4.*t9.*u6);t2.*(-t3.*t4.*t15+t9.*t32.*u2.*u4+t4.*t41.*u4.*u6-t4.*t41.*u2.*u16.*(2.0./3.0)+param1.*t4.*t9.*t63.*t64.*(t32.*u1.*(u1.*(t46-t91)+t9.*u4.*u6)-t3.*t32.*u4.*u6));zero;zero;zero;zero;-t2.*(u2.*(t3+t3.*t32)-param1.*t4.*t9.*t63.*u6);zero;t2.*t4.*t9.*u2.*2.0;t2.*t4.*t9.*u3.*2.0;t94;t2.*(t4.*t8.*t41.*(4.0./3.0)+t4.*t10.*t41+t4.*t37.*t41+param1.*t4.*t9.*t63.*t64.*(t32.*t40+t32.*t76.*u1));zero;t2.*t3.*t4.*-2.0;zero;zero;-t2.*(t4.*t9.*u2.*(4.0./3.0)-param1.*t4.*t9.*t63.*u2);zero;zero;t2.*t3.*t4.*-2.0;zero;-t2.*(t4.*t9.*u3-param1.*t4.*t9.*t63.*u3);zero;zero;zero;-t2.*t3.*t4;-t2.*(t4.*t9.*u4-param1.*t4.*t9.*t63.*u4);zero;zero;zero;zero;-param1.*t2.*t3.*t4.*t63;zero;t4.*t9.*t93.*u3.*-2.0;t4.*t9.*t93.*u2.*2.0;zero;t4.*t41.*t93.*u2.*u3.*(1.0./3.0);zero;zero;t3.*t4.*t93.*-2.0;zero;-t4.*t9.*t93.*u3;zero;t3.*t4.*t93.*2.0;zero;zero;t4.*t9.*t93.*u2.*(2.0./3.0);zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;t2.*t4.*t9.*u2;t2.*t4.*t41.*u2.*u4.*(1.0./3.0);zero;zero;zero;-t2.*t3.*t4;-t94;zero;zero;zero;zero;zero;zero;zero;zero;zero;t2.*t4.*t9.*u2.*(2.0./3.0);zero;zero;zero;zero;zero];
end
s = reshape(s,ng,nch);
s_udg = reshape(s_udg,ng,nch,nc);
