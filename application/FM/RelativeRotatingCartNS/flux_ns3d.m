function [f,f_udg] = flux_ns3d(pg,udg,param,time)
%FLUX_NS3D
%    [F,F_UDG] = FLUX_NS3D(PG,UDG,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    16-Dec-2019 16:05:09
[ng,nc] = size(udg);
nch = 5;
nd = 3;
one = ones(ng,1);
param1 = param{1};
param2 = param{2};
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
u14 = udg(:,14);
u15 = udg(:,15);
u16 = udg(:,16);
u17 = udg(:,17);
u18 = udg(:,18);
u19 = udg(:,19);
u20 = udg(:,20);
x1 = pg(:,1);
x2 = pg(:,2);
zero = zeros(ng,1);
t2 = 1.0./u1;
t3 = param2.^2;
t4 = u2.^2;
t5 = 1.0./u1.^2;
t6 = 1.0./param3;
t7 = param1-1.0;
t8 = x1.^2;
t9 = t3.*t8.*(1.0./2.0);
t10 = x2.^2;
t11 = t3.*t10.*(1.0./2.0);
t12 = t9+t11;
t13 = t12.*u1;
t14 = t4.*t5.*(1.0./2.0);
t15 = u3.^2;
t16 = t5.*t15.*(1.0./2.0);
t17 = u4.^2;
t18 = t5.*t17.*(1.0./2.0);
t19 = t14+t16+t18;
t37 = t19.*u1;
t20 = t13-t37+u5;
t39 = t2.*u3.*u6;
t21 = -t39+u8;
t22 = t2.*t21;
t41 = t2.*u2.*u11;
t23 = -t41+u12;
t24 = t2.*t23;
t25 = t22+t24;
t40 = t2.*u4.*u6;
t26 = -t40+u9;
t27 = t2.*t26;
t61 = t2.*u2.*u16;
t28 = -t61+u17;
t29 = t2.*t28;
t30 = t27+t29;
t38 = t2.*u2.*u6;
t31 = -t38+u7;
t45 = t2.*u3.*u11;
t32 = -t45+u13;
t33 = t2.*t32;
t46 = t2.*u4.*u16;
t34 = -t46+u19;
t35 = t2.*t34;
t77 = t2.*t31.*2.0;
t36 = t33+t35-t77;
t42 = t6.*t25;
t43 = t2.*u2.*u3;
t44 = t42+t43;
t47 = t7.*t20;
t48 = t2.*u5;
t49 = t2.*t7.*t20;
t50 = t48+t49;
t60 = t2.*u4.*u11;
t51 = -t60+u14;
t52 = t2.*t51;
t65 = t2.*u3.*u16;
t53 = -t65+u18;
t54 = t2.*t53;
t55 = t52+t54;
t56 = t2.*t31;
t109 = t2.*t32.*2.0;
t57 = t35+t56-t109;
t58 = 1.0./param4;
t59 = 1.0./t7;
t62 = t6.*t30;
t63 = t2.*u2.*u4;
t64 = t62+t63;
t66 = t6.*t55;
t67 = t2.*u3.*u4;
t68 = t66+t67;
t137 = t2.*t34.*2.0;
t69 = t33+t56-t137;
t70 = 1.0./u1.^3;
t71 = t4.*t70;
t72 = t15.*t70;
t73 = t17.*t70;
t74 = t71+t72+t73;
t75 = t74.*u1;
t76 = t9+t11-t14-t16-t18+t75;
t78 = t5.*t21;
t79 = t5.*t23;
t101 = t70.*u3.*u6;
t102 = t70.*u2.*u11;
t80 = t78+t79-t101-t102;
t81 = t5.*t26;
t82 = t5.*t28;
t128 = t70.*u4.*u6;
t129 = t70.*u2.*u16;
t83 = t81+t82-t128-t129;
t84 = t5.*t31.*2.0;
t85 = t70.*u3.*u11;
t86 = t70.*u4.*u16;
t105 = t5.*t34;
t135 = t5.*t32;
t87 = t84+t85+t86-t105-t135-t70.*u2.*u6.*2.0;
t88 = t7.*t20.*u6;
t89 = t5.*t31.*u2;
t90 = t5.*t21.*u3;
t91 = t5.*t26.*u4;
t92 = t89+t90+t91;
t93 = t92.*u1;
t94 = t19.*u6;
t95 = t3.*u1.*x1;
t99 = t12.*u6;
t96 = t93+t94+t95-t99-u10;
t97 = t7.*t96.*u1;
t98 = t88+t97;
t100 = 1.0./u1.^4;
t103 = -t6.*t80-t5.*u2.*u3;
t104 = t7.*t76;
t106 = t5.*u5;
t107 = t5.*t7.*t20;
t136 = t2.*t7.*t76;
t108 = t106+t107-t136;
t110 = t5.*t51;
t111 = t5.*t53;
t131 = t70.*u4.*u11;
t132 = t70.*u3.*u16;
t112 = t110+t111-t131-t132;
t113 = t5.*t32.*2.0;
t114 = t70.*u2.*u6;
t134 = t5.*t31;
t115 = t86-t105+t113+t114-t134-t70.*u3.*u11.*2.0;
t116 = t7.*t20.*u11;
t117 = t5.*t23.*u2;
t118 = t5.*t32.*u3;
t119 = t5.*t51.*u4;
t120 = t117+t118+t119;
t121 = t120.*u1;
t122 = t19.*u11;
t123 = t3.*u1.*x2;
t127 = t12.*u11;
t124 = t121+t122+t123-t127-u15;
t125 = t7.*t124.*u1;
t126 = t116+t125;
t130 = -t6.*t83-t5.*u2.*u4;
t133 = -t6.*t112-t5.*u3.*u4;
t138 = t5.*t34.*2.0;
t139 = t85+t114-t134-t135+t138-t70.*u4.*u16.*2.0;
t140 = t5.*t28.*u2;
t141 = t5.*t53.*u3;
t142 = t5.*t34.*u4;
t143 = t140+t141+t142;
t144 = t12.*u16;
t147 = t143.*u1;
t148 = t19.*u16;
t145 = t144-t147-t148+u20;
t146 = t7.*t20.*u16;
t149 = t146-t7.*t145.*u1;
f = [u2;t47+t2.*t4-t6.*t36.*(2.0./3.0);t44;t64;t50.*u2+t2.*t6.*t25.*u3+t2.*t6.*t30.*u4-t2.*t6.*t36.*u2.*(2.0./3.0)-param1.*t5.*t6.*t58.*t59.*t98;u3;t44;t47+t2.*t15-t6.*t57.*(2.0./3.0);t68;t50.*u3+t2.*t6.*t25.*u2+t2.*t6.*t55.*u4-t2.*t6.*t57.*u3.*(2.0./3.0)-param1.*t5.*t6.*t58.*t59.*t126;u4;t64;t68;t47+t2.*t17-t6.*t69.*(2.0./3.0);t50.*u4+t2.*t6.*t30.*u2+t2.*t6.*t55.*u3-t2.*t6.*t69.*u4.*(2.0./3.0)-param1.*t5.*t6.*t58.*t59.*t149];
if nargout > 1
    t150 = t2.*u3;
    t166 = t5.*t6.*u11;
    t151 = t150-t166;
    t152 = t2.*u4;
    t159 = t5.*t6.*u16;
    t153 = t152-t159;
    t154 = t5.*t6.*u6.*(2.0./3.0);
    t155 = t154-t2.*t7.*u2;
    t156 = t2.*t6.*t25;
    t157 = t2.*u2;
    t162 = t5.*t6.*u6;
    t158 = t157-t162;
    t160 = t5.*t6.*u11.*(2.0./3.0);
    t161 = t160-t2.*t7.*u3;
    t163 = t2.*t6.*t30;
    t164 = t5.*t6.*u16.*(2.0./3.0);
    t168 = t2.*t7.*u4;
    t165 = t164-t168;
    t167 = t2.*t6.*t55;
    t169 = one.*t7;
    t170 = t2.*t7;
    t171 = t2+t170;
    t172 = t5.*t6.*u2.*(2.0./3.0);
    t173 = t2.*t6;
    t174 = t5.*t6.*u2;
    t175 = t5.*t6.*u4;
    t183 = t7.*t76.*u1;
    t176 = t47-t183;
    t177 = t5.*t6.*u3.*(2.0./3.0);
    t178 = t5.*t6.*u3;
    t179 = t2.*t6.*(4.0./3.0);
    t186 = param1.*t5.*t6.*t58.*u4;
    t180 = t175-t186;
    t181 = param1.*t2.*t6.*t58;
    t182 = t5.*t6.*u4.*(2.0./3.0);
    t184 = t174-param1.*t5.*t6.*t58.*u2;
    t185 = t178-param1.*t5.*t6.*t58.*u3;
    f_udg = [zero;t104-t4.*t5-t6.*t87.*(2.0./3.0);t103;t130;-t108.*u2-t5.*t6.*t25.*u3-t5.*t6.*t30.*u4+t5.*t6.*t36.*u2.*(2.0./3.0)-t2.*t6.*t80.*u3-t2.*t6.*t83.*u4-t2.*t6.*t87.*u2.*(2.0./3.0)+param1.*t6.*t58.*t59.*t70.*t98.*2.0-param1.*t5.*t6.*t58.*t59.*(t7.*t96+t7.*u1.*(t89+t90+t91-t74.*u6+t3.*x1-u1.*(t21.*t70.*u3.*2.0+t26.*t70.*u4.*2.0+t31.*t70.*u2.*2.0-t4.*t100.*u6-t15.*t100.*u6-t17.*t100.*u6))+t7.*t76.*u6);zero;t103;t104-t5.*t15-t6.*t115.*(2.0./3.0);t133;-t108.*u3-t5.*t6.*t25.*u2-t5.*t6.*t55.*u4+t5.*t6.*t57.*u3.*(2.0./3.0)-t2.*t6.*t80.*u2-t2.*t6.*t112.*u4-t2.*t6.*t115.*u3.*(2.0./3.0)+param1.*t6.*t58.*t59.*t70.*t126.*2.0-param1.*t5.*t6.*t58.*t59.*(t7.*t124+t7.*u1.*(t117+t118+t119-t74.*u11+t3.*x2-u1.*(t23.*t70.*u2.*2.0+t32.*t70.*u3.*2.0-t4.*t100.*u11+t51.*t70.*u4.*2.0-t15.*t100.*u11-t17.*t100.*u11))+t7.*t76.*u11);zero;t130;t133;t104-t5.*t17-t6.*t139.*(2.0./3.0);-t108.*u4-t5.*t6.*t30.*u2-t5.*t6.*t55.*u3+t5.*t6.*t69.*u4.*(2.0./3.0)-t2.*t6.*t83.*u2-t2.*t6.*t112.*u3-t2.*t6.*t139.*u4.*(2.0./3.0)+param1.*t6.*t58.*t59.*t70.*t149.*2.0-param1.*t5.*t6.*t58.*t59.*(-t7.*t145+t7.*t76.*u16+t7.*u1.*(t140+t141+t142-t74.*u16-u1.*(t28.*t70.*u2.*2.0+t34.*t70.*u4.*2.0-t4.*t100.*u16+t53.*t70.*u3.*2.0-t15.*t100.*u16-t17.*t100.*u16)));one;t2.*u2.*2.0-t2.*t7.*u2-t5.*t6.*u6.*(4.0./3.0);t151;t153;t48+t49-t4.*t5.*t7-t2.*t6.*t36.*(2.0./3.0)-t6.*t70.*u2.*u6.*(4.0./3.0)-t6.*t70.*u3.*u11-t6.*t70.*u4.*u16+param1.*t5.*t6.*t58.*t59.*(t7.*u1.*(u1.*(t114-t134)-t5.*u2.*u6)+t2.*t7.*u2.*u6);zero;t151;t155;zero;t156-t5.*t7.*u2.*u3+t6.*t70.*u3.*u6.*(2.0./3.0)-t6.*t70.*u2.*u11-param1.*t5.*t6.*t58.*t59.*(t7.*u1.*(u1.*(t79-t102)+t5.*u2.*u11)-t2.*t7.*u2.*u11);zero;t153;zero;t155;t163-t5.*t7.*u2.*u4+t6.*t70.*u4.*u6.*(2.0./3.0)-t6.*t70.*u2.*u16-param1.*t5.*t6.*t58.*t59.*(t7.*u1.*(u1.*(t82-t129)+t5.*u2.*u16)-t2.*t7.*u2.*u16);zero;t161;t158;zero;t156-t5.*t7.*u2.*u3-t6.*t70.*u3.*u6+t6.*t70.*u2.*u11.*(2.0./3.0)-param1.*t5.*t6.*t58.*t59.*(t7.*u1.*(u1.*(t78-t101)+t5.*u3.*u6)-t2.*t7.*u3.*u6);one;t158;t2.*u3.*2.0-t2.*t7.*u3-t5.*t6.*u11.*(4.0./3.0);t153;t48+t49-t5.*t7.*t15-t2.*t6.*t57.*(2.0./3.0)-t6.*t70.*u2.*u6-t6.*t70.*u3.*u11.*(4.0./3.0)-t6.*t70.*u4.*u16+param1.*t5.*t6.*t58.*t59.*(t7.*u1.*(u1.*(t85-t135)-t5.*u3.*u11)+t2.*t7.*u3.*u11);zero;zero;t153;t161;t167-t5.*t7.*u3.*u4+t6.*t70.*u4.*u11.*(2.0./3.0)-t6.*t70.*u3.*u16-param1.*t5.*t6.*t58.*t59.*(t7.*u1.*(u1.*(t111-t132)+t5.*u3.*u16)-t2.*t7.*u3.*u16);zero;t165;zero;t158;t163-t5.*t7.*u2.*u4-t6.*t70.*u4.*u6+t6.*t70.*u2.*u16.*(2.0./3.0)-param1.*t5.*t6.*t58.*t59.*(t7.*u1.*(u1.*(t81-t128)+t5.*u4.*u6)-t2.*t7.*u4.*u6);zero;zero;t165;t151;t167-t5.*t7.*u3.*u4-t6.*t70.*u4.*u11+t6.*t70.*u3.*u16.*(2.0./3.0)-param1.*t5.*t6.*t58.*t59.*(t7.*u1.*(u1.*(t110-t131)+t5.*u4.*u11)-t2.*t7.*u4.*u11);one;t158;t151;-t168+t2.*u4.*2.0-t5.*t6.*u16.*(4.0./3.0);t48+t49-t5.*t7.*t17-t2.*t6.*t69.*(2.0./3.0)-t6.*t70.*u2.*u6-t6.*t70.*u3.*u11-t6.*t70.*u4.*u16.*(4.0./3.0)+param1.*t5.*t6.*t58.*t59.*(t7.*u1.*(u1.*(t86-t105)-t5.*u4.*u16)+t2.*t7.*u4.*u16);zero;t169;zero;zero;t171.*u2-param1.*t5.*t6.*t58.*u6;zero;zero;t169;zero;t171.*u3-param1.*t5.*t6.*t58.*u11;zero;zero;zero;t169;t171.*u4-param1.*t5.*t6.*t58.*u16;zero;t5.*t6.*u2.*(-4.0./3.0);-t5.*t6.*u3;-t5.*t6.*u4;t4.*t6.*t70.*(-4.0./3.0)-t6.*t15.*t70-t6.*t17.*t70-param1.*t5.*t6.*t58.*t59.*t176;zero;-t5.*t6.*u3;t172;zero;t6.*t70.*u2.*u3.*(-1.0./3.0);zero;-t5.*t6.*u4;zero;t172;t6.*t70.*u2.*u4.*(-1.0./3.0);zero;t179;zero;zero;t5.*t6.*u2.*(4.0./3.0)-param1.*t5.*t6.*t58.*u2;zero;zero;t2.*t6.*(-2.0./3.0);zero;t5.*t6.*u3.*(-2.0./3.0);zero;zero;zero;t2.*t6.*(-2.0./3.0);t5.*t6.*u4.*(-2.0./3.0);zero;zero;t173;zero;t185;zero;t173;zero;zero;t174;zero;zero;zero;zero;zero;zero;zero;zero;t173;t180;zero;zero;zero;zero;zero;zero;t173;zero;zero;t174;zero;zero;zero;zero;t181;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;t177;-t174;zero;t6.*t70.*u2.*u3.*(-1.0./3.0);zero;-t174;t5.*t6.*u3.*(-4.0./3.0);-t175;-t4.*t6.*t70-t6.*t15.*t70.*(4.0./3.0)-t6.*t17.*t70-param1.*t5.*t6.*t58.*t59.*t176;zero;zero;-t175;t177;t6.*t70.*u3.*u4.*(-1.0./3.0);zero;zero;t173;zero;t178;zero;t173;zero;zero;t184;zero;zero;zero;zero;zero;zero;t2.*t6.*(-2.0./3.0);zero;zero;-t172;zero;zero;t179;zero;t5.*t6.*u3.*(4.0./3.0)-param1.*t5.*t6.*t58.*u3;zero;zero;zero;t2.*t6.*(-2.0./3.0);t5.*t6.*u4.*(-2.0./3.0);zero;zero;zero;zero;zero;zero;zero;zero;t173;t180;zero;zero;t173;zero;t178;zero;zero;zero;zero;zero;zero;zero;zero;zero;t181;zero;zero;zero;zero;zero;zero;t182;zero;-t174;t6.*t70.*u2.*u4.*(-1.0./3.0);zero;zero;t182;-t178;t6.*t70.*u3.*u4.*(-1.0./3.0);zero;-t174;-t178;t5.*t6.*u4.*(-4.0./3.0);-t4.*t6.*t70-t6.*t15.*t70-t6.*t17.*t70.*(4.0./3.0)-param1.*t5.*t6.*t58.*t59.*t176;zero;zero;zero;t173;t175;zero;zero;zero;zero;zero;zero;t173;zero;zero;t184;zero;zero;zero;zero;zero;zero;zero;zero;t173;t175;zero;zero;t173;zero;t185;zero;t2.*t6.*(-2.0./3.0);zero;zero;-t172;zero;zero;t2.*t6.*(-2.0./3.0);zero;-t177;zero;zero;zero;t179;-t186+t5.*t6.*u4.*(4.0./3.0);zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;t181];
end
f = reshape(f,ng,nch,nd);
f_udg = reshape(f_udg,ng,nch,nd,nc);
