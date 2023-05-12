function [fh,fh_udg,fh_uh] = fhat_ns3d(nl,pg,udg,uh,param,time)
%FHAT_NS3D
%    [FH,FH_UDG,FH_UH] = FHAT_NS3D(NL,PG,UDG,UH,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    16-Dec-2019 16:05:15
[ng,nc] = size(udg);
nch = 5;
nd = 3;
nl1 = nl(:,1);
nl2 = nl(:,2);
nl3 = nl(:,3);
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
u13 = udg(:,13);
u14 = udg(:,14);
u15 = udg(:,15);
u16 = udg(:,16);
u17 = udg(:,17);
u18 = udg(:,18);
u19 = udg(:,19);
u20 = udg(:,20);
uh5 = uh(:,5);
x1 = pg(:,1);
x2 = pg(:,2);
zero = zeros(ng,1);
t2 = param1-1.0;
t3 = param2.^2;
t4 = 1.0./u1.^2;
t5 = 1.0./u1;
t6 = x1.^2;
t7 = t3.*t6.*(1.0./2.0);
t8 = x2.^2;
t9 = t3.*t8.*(1.0./2.0);
t10 = t7+t9;
t11 = u2.^2;
t12 = t4.*t11.*(1.0./2.0);
t13 = u3.^2;
t14 = t4.*t13.*(1.0./2.0);
t15 = u4.^2;
t16 = t4.*t15.*(1.0./2.0);
t17 = t12+t14+t16;
t18 = param5.^2;
t19 = t10.*u1;
t21 = t17.*u1;
t20 = t19-t21+u5;
t22 = 1.0./param3;
t23 = 1.0./param4;
t24 = 1.0./param5.^2;
t25 = 1.0./t2;
t39 = t5.*u2.*u6;
t26 = -t39+u7;
t27 = t4.*t26.*u2;
t40 = t5.*u3.*u6;
t28 = -t40+u8;
t29 = t4.*t28.*u3;
t41 = t5.*u4.*u6;
t30 = -t41+u9;
t31 = t4.*t30.*u4;
t32 = t27+t29+t31;
t33 = t32.*u1;
t34 = t17.*u6;
t35 = t3.*u1.*x1;
t73 = t10.*u6;
t36 = t33+t34+t35-t73-u10;
t37 = 1.0./u1.^3;
t38 = 1.0./u1.^4;
t42 = t11.*t37;
t43 = t13.*t37;
t44 = t15.*t37;
t45 = t42+t43+t44;
t57 = t5.*u2.*u11;
t46 = -t57+u12;
t47 = t4.*t46.*u2;
t58 = t5.*u3.*u11;
t48 = -t58+u13;
t49 = t4.*t48.*u3;
t59 = t5.*u4.*u11;
t50 = -t59+u14;
t51 = t4.*t50.*u4;
t52 = t47+t49+t51;
t53 = t52.*u1;
t54 = t17.*u11;
t55 = t3.*u1.*x2;
t77 = t10.*u11;
t56 = t53+t54+t55-t77-u15;
t60 = t45.*u1;
t61 = t7+t9-t12-t14-t16+t60;
t80 = t5.*u2.*u16;
t62 = -t80+u17;
t63 = t4.*t62.*u2;
t81 = t5.*u3.*u16;
t64 = -t81+u18;
t65 = t4.*t64.*u3;
t82 = t5.*u4.*u16;
t66 = -t82+u19;
t67 = t4.*t66.*u4;
t68 = t63+t65+t67;
t69 = t10.*u16;
t83 = t68.*u1;
t84 = t17.*u16;
t70 = t69-t83-t84+u20;
t71 = t2.*t20.*u16-t2.*t70.*u1;
t72 = t2.*t20.*u6;
t74 = t2.*t36.*u1;
t75 = t72+t74;
t76 = t2.*t20.*u11;
t78 = t2.*t56.*u1;
t79 = t76+t78;
fh = param6.*(u5-uh5)-t22.*t23.*t24.*t25.*(nl3.*param1.*t4.*t18.*t71+nl1.*param1.*t4.*t18.*t75+nl2.*param1.*t4.*t18.*t79);
if nargout > 1
    t85 = t2.*t20;
    t87 = t2.*t61.*u1;
    t86 = t85-t87;
    fh_udg = [-t22.*t23.*t24.*t25.*(nl3.*param1.*t4.*t18.*(-t2.*t70+t2.*t61.*u16+t2.*u1.*(t63+t65+t67-t45.*u16+u1.*(t11.*t38.*u16+t13.*t38.*u16+t15.*t38.*u16-t37.*t62.*u2.*2.0-t37.*t64.*u3.*2.0-t37.*t66.*u4.*2.0)))-nl3.*param1.*t18.*t37.*t71.*2.0-nl1.*param1.*t18.*t37.*t75.*2.0-nl2.*param1.*t18.*t37.*t79.*2.0+nl1.*param1.*t4.*t18.*(t2.*t36+t2.*u1.*(t27+t29+t31-t45.*u6+t3.*x1+u1.*(t11.*t38.*u6+t13.*t38.*u6+t15.*t38.*u6-t26.*t37.*u2.*2.0-t28.*t37.*u3.*2.0-t30.*t37.*u4.*2.0))+t2.*t61.*u6)+nl2.*param1.*t4.*t18.*(t2.*t56+t2.*u1.*(t47+t49+t51-t45.*u11+t3.*x2+u1.*(t11.*t38.*u11+t13.*t38.*u11+t15.*t38.*u11-t37.*t46.*u2.*2.0-t37.*t48.*u3.*2.0-t37.*t50.*u4.*2.0))+t2.*t61.*u11));-t22.*t23.*t24.*t25.*(nl1.*param1.*t4.*t18.*(t2.*u1.*(u1.*(t4.*t26-t37.*u2.*u6)+t4.*u2.*u6)-t2.*t5.*u2.*u6)+nl2.*param1.*t4.*t18.*(t2.*u1.*(u1.*(t4.*t46-t37.*u2.*u11)+t4.*u2.*u11)-t2.*t5.*u2.*u11)+nl3.*param1.*t4.*t18.*(t2.*u1.*(u1.*(t4.*t62-t37.*u2.*u16)+t4.*u2.*u16)-t2.*t5.*u2.*u16));-t22.*t23.*t24.*t25.*(nl1.*param1.*t4.*t18.*(t2.*u1.*(u1.*(t4.*t28-t37.*u3.*u6)+t4.*u3.*u6)-t2.*t5.*u3.*u6)+nl2.*param1.*t4.*t18.*(t2.*u1.*(u1.*(t4.*t48-t37.*u3.*u11)+t4.*u3.*u11)-t2.*t5.*u3.*u11)+nl3.*param1.*t4.*t18.*(t2.*u1.*(u1.*(t4.*t64-t37.*u3.*u16)+t4.*u3.*u16)-t2.*t5.*u3.*u16));-t22.*t23.*t24.*t25.*(nl1.*param1.*t4.*t18.*(t2.*u1.*(u1.*(t4.*t30-t37.*u4.*u6)+t4.*u4.*u6)-t2.*t5.*u4.*u6)+nl2.*param1.*t4.*t18.*(t2.*u1.*(u1.*(t4.*t50-t37.*u4.*u11)+t4.*u4.*u11)-t2.*t5.*u4.*u11)+nl3.*param1.*t4.*t18.*(t2.*u1.*(u1.*(t4.*t66-t37.*u4.*u16)+t4.*u4.*u16)-t2.*t5.*u4.*u16));param6-t22.*t23.*t24.*t25.*(nl1.*param1.*t2.*t4.*t18.*u6+nl2.*param1.*t2.*t4.*t18.*u11+nl3.*param1.*t2.*t4.*t18.*u16);-nl1.*param1.*t4.*t22.*t23.*t25.*t86;-nl1.*param1.*t4.*t22.*t23.*u2;-nl1.*param1.*t4.*t22.*t23.*u3;-nl1.*param1.*t4.*t22.*t23.*u4;nl1.*param1.*t5.*t22.*t23;-nl2.*param1.*t4.*t22.*t23.*t25.*t86;-nl2.*param1.*t4.*t22.*t23.*u2;-nl2.*param1.*t4.*t22.*t23.*u3;-nl2.*param1.*t4.*t22.*t23.*u4;nl2.*param1.*t5.*t22.*t23;-nl3.*param1.*t4.*t22.*t23.*t25.*t86;-nl3.*param1.*t4.*t22.*t23.*u2;-nl3.*param1.*t4.*t22.*t23.*u3;-nl3.*param1.*t4.*t22.*t23.*u4;nl3.*param1.*t5.*t22.*t23];
end
if nargout > 2
    fh_uh = [zero;zero;zero;zero;-one.*param6];
end
fh = reshape(fh,ng,1);
fh_udg = reshape(fh_udg,ng,1,nc);
fh_uh = reshape(fh_uh,ng,1,nch);

