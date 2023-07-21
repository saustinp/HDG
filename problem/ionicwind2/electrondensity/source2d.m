function [s,s_udg] = source2d(p,udg,param,time)
%SOURCE2D
%    [S,S_UDG] = SOURCE2D(P,UDG,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    21-Jul-2023 01:06:56
[ng,nc] = size(udg);
nch = 4;
nd = 2;
one = ones(ng,1);
param5 = param{5};
param6 = param{6};
param11 = param{11};
param12 = param{12};
param15 = param{15};
param16 = param{16};
param18 = param{18};
param19 = param{19};
u1 = udg(:,1);
u2 = udg(:,2);
u3 = udg(:,3);
u8 = udg(:,8);
u12 = udg(:,12);
x1 = p(:,1);
x3 = p(:,3);
x4 = p(:,4);
zero = zeros(ng,1);
t2 = param5.*u1;
t3 = param6.*u2;
t4 = one.*x1;
t5 = u8+x3;
t6 = u12+x4;
t7 = log(1.0e+1);
t12 = 1.0./param11;
t13 = 1.0./param18;
t14 = 1.0./param19;
t15 = -t4;
t16 = t5.^2;
t17 = t6.^2;
t18 = 1.0./t7;
t24 = t2+t3;
t26 = param5.*param12.*t12.*t14.*u3;
t27 = param6.*param12.*t12.*t14.*u3.*x1;
t19 = t18.^2;
t20 = t18.^3;
t25 = t16+t17;
t30 = -t27;
t21 = t19.^2;
t28 = 1.0./t25;
t29 = sqrt(t25);
t31 = 1.0./t29;
t32 = param15.*t13.*t29.*1.0e+21;
t34 = param15.*t13.*t29.*1.0e+24;
t43 = t5.*t18.*t28.*3.042512694247558e+2;
t44 = t6.*t18.*t28.*3.042512694247558e+2;
t46 = t5.*t18.*t28.*1.683358794346262e+2;
t47 = t6.*t18.*t28.*1.683358794346262e+2;
t33 = log(t32);
t38 = t34-2.009e+4;
t35 = t33.^2;
t36 = t33.^3;
t39 = tanh(t38);
t40 = t18.*t33.*1.0e+3;
t56 = t18.*t33.*1.044327759364851;
t57 = t18.*t33.*1.683358794346262e+2;
t58 = t18.*t33.*3.042512694247558e+2;
t65 = t18.*t33.*3.4808749e-1;
t72 = t5.*t19.*t28.*t33.*3.725749808547279e+2;
t73 = t6.*t19.*t28.*t33.*3.725749808547279e+2;
t74 = t5.*t19.*t28.*t33.*2.072993272522841e+2;
t75 = t6.*t19.*t28.*t33.*2.072993272522841e+2;
t37 = t35.^2;
t41 = t39.^2;
t42 = t39./2.0;
t48 = t40-1.1e+3;
t59 = -t57;
t60 = t19.*t35.*1.862874904273639e+2;
t62 = t20.*t36.*5.150075493521921e+1;
t63 = t20.*t36.*2.845690675675998e+1;
t64 = t19.*t35.*1.036496636261421e+2;
t69 = -t65;
t76 = -t72;
t77 = -t73;
t78 = -t74;
t79 = -t75;
t80 = t5.*t21.*t28.*t36.*2.144751826853317e+1;
t81 = t6.*t21.*t28.*t36.*2.144751826853317e+1;
t82 = t5.*t21.*t28.*t36.*1.177095445972644e+1;
t83 = t6.*t21.*t28.*t36.*1.177095445972644e+1;
t84 = t5.*t20.*t28.*t35.*8.537072027027995e+1;
t85 = t6.*t20.*t28.*t35.*8.537072027027995e+1;
t86 = t5.*t20.*t28.*t35.*1.545022648056576e+2;
t87 = t6.*t20.*t28.*t35.*1.545022648056576e+2;
t92 = t56+4.012938813577186e+1;
t45 = t41-1.0;
t49 = tanh(t48);
t50 = t42+1.0./2.0;
t61 = -t60;
t66 = t21.*t37.*2.94273861493161;
t67 = -t63;
t68 = t21.*t37.*5.361879567133292;
t71 = 1.0e+1.^t69;
t88 = -t80;
t89 = -t81;
t90 = -t82;
t91 = -t83;
t51 = t49.^2;
t52 = t49./2.0;
t70 = -t68;
t101 = t59+t64+t66+t67+1.251804446871849e+2;
t103 = t43+t76+t86+t88;
t104 = t44+t77+t87+t89;
t105 = t46+t78+t84+t90;
t106 = t47+t79+t85+t91;
t53 = t51-1.0;
t54 = t52+1.0./2.0;
t55 = t52-1.0./2.0;
t100 = t58+t61+t62+t70-2.123274758440036e+2;
t93 = t5.*t18.*t28.*t55.*1.044327759364851;
t94 = t6.*t18.*t28.*t55.*1.044327759364851;
t95 = t55.*t92;
t96 = t5.*t18.*t28.*t53.*t92.*5.0e+2;
t97 = t6.*t18.*t28.*t53.*t92.*5.0e+2;
t102 = 1.0e+1.^t100;
t109 = t54.*t101;
t111 = t54.*t105;
t112 = t54.*t106;
t113 = t5.*t18.*t28.*t53.*t101.*5.0e+2;
t114 = t6.*t18.*t28.*t53.*t101.*5.0e+2;
t98 = -t96;
t99 = -t97;
t107 = param18.*t50.*t102;
t110 = -t109;
t108 = -t107;
t115 = t95+t110;
t119 = t93+t98+t111+t113;
t120 = t94+t99+t112+t114;
t116 = 1.0e+1.^t115;
t117 = param18.*t116;
s = [-x1.*(param12.*t2.*t12.*t14.*u3-param16.*t13.*t14.*t29.*t71.*u1.*(t107-t117).*5.849682301203802e+24);-x1.*(param12.*t3.*t12.*t14.*u3-param16.*t14.*t29.*t71.*t116.*u1.*5.849682301203802e+24);-x1.*(param12.*t12.*t14.*t24.*u3-param16.*t14.*t29.*t50.*t71.*t102.*u1.*5.849682301203802e+24);-x1.*(u1+u2-u3)];
if nargout > 1
    mt1 = [-x1.*(t26-param16.*t13.*t14.*t29.*t71.*(t107-t117).*5.849682301203802e+24);param16.*t14.*t29.*t71.*t116.*x1.*5.849682301203802e+24;-x1.*(t26-param16.*t14.*t29.*t50.*t71.*t102.*5.849682301203802e+24);t15;zero;t30;t30;t15;-param12.*t2.*t12.*t14.*x1;-param12.*t3.*t12.*t14.*x1;-param12.*t12.*t14.*t24.*x1;t4;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero];
    mt2 = [-x1.*(param16.*t13.*t14.*t29.*t71.*u1.*(t7.*t103.*t108+t7.*t117.*t119+param15.*t5.*t31.*t45.*t102.*5.0e+23).*5.849682301203802e+24-param16.*t5.*t13.*t14.*t31.*t71.*u1.*(t107-t117).*3.813481071680347e+24)];
    mt3 = [x1.*(param16.*t5.*t14.*t31.*t71.*t116.*u1.*3.813481071680347e+24+param16.*t7.*t14.*t29.*t71.*t116.*t119.*u1.*5.849682301203802e+24)];
    mt4 = [x1.*(param16.*t5.*t14.*t31.*t50.*t71.*t102.*u1.*3.813481071680347e+24-param15.*param16.*t5.*t13.*t14.*t45.*t71.*t102.*u1.*2.924841150601901e+48+param16.*t7.*t14.*t29.*t50.*t71.*t102.*t103.*u1.*5.849682301203802e+24);zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero;zero];
    mt5 = [-x1.*(param16.*t13.*t14.*t29.*t71.*u1.*(t7.*t104.*t108+t7.*t117.*t120+param15.*t6.*t31.*t45.*t102.*5.0e+23).*5.849682301203802e+24-param16.*t6.*t13.*t14.*t31.*t71.*u1.*(t107-t117).*3.813481071680347e+24)];
    mt6 = [x1.*(param16.*t6.*t14.*t31.*t71.*t116.*u1.*3.813481071680347e+24+param16.*t7.*t14.*t29.*t71.*t116.*t120.*u1.*5.849682301203802e+24)];
    mt7 = [x1.*(param16.*t6.*t14.*t31.*t50.*t71.*t102.*u1.*3.813481071680347e+24-param15.*param16.*t6.*t13.*t14.*t45.*t71.*t102.*u1.*2.924841150601901e+48+param16.*t7.*t14.*t29.*t50.*t71.*t102.*t104.*u1.*5.849682301203802e+24);zero];
    s_udg = [mt1;mt2;mt3;mt4;mt5;mt6;mt7];
end
s = reshape(s,ng,nch);
s_udg = reshape(s_udg,ng,nch,nc);
end