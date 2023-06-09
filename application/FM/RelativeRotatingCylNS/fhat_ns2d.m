function [fh,fh_udg,fh_uh] = fhat_ns2d(nl,pg,udg,uh,param,time)
%FHAT_NS2D
%    [FH,FH_UDG,FH_UH] = FHAT_NS2D(NL,PG,UDG,UH,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    14-Dec-2019 22:17:52
[ng,nc] = size(udg);
nch = 4;
nd = 2;
nl1 = nl(:,1);
nl2 = nl(:,2);
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
uh4 = uh(:,4);
x1 = pg(:,1);
zero = zeros(ng,1);
t2 = param1-1.0;
t3 = 1.0./u1.^2;
t4 = 1.0./u1;
t5 = param2.^2;
t6 = u2.^2;
t7 = t3.*t6.*(1.0./2.0);
t8 = u3.^2;
t9 = t3.*t8.*(1.0./2.0);
t10 = t7+t9;
t11 = x1.^2;
t12 = param5.^2;
t13 = t5.*t11.*u1.*(1.0./2.0);
t29 = t10.*u1;
t14 = t13-t29+u4;
t15 = 1.0./param3;
t16 = 1.0./param4;
t17 = 1.0./param5.^2;
t18 = 1.0./t2;
t19 = t10.*u5;
t32 = t4.*u2.*u5;
t20 = -t32+u6;
t21 = t3.*t20.*u2;
t33 = t4.*u3.*u5;
t22 = -t33+u7;
t23 = t3.*t22.*u3;
t24 = t21+t23;
t25 = t24.*u1;
t26 = t5.*u1.*x1;
t34 = t5.*t11.*u5.*(1.0./2.0);
t27 = t19+t25+t26-t34-u8;
t28 = t2.*t27.*u1;
t30 = t2.*t14.*u5;
t31 = t28+t30;
t35 = 1.0./u1.^3;
t36 = t6.*t35;
t37 = t8.*t35;
t38 = t36+t37;
t39 = 1.0./u1.^4;
t40 = 1.0./x1.^2;
t49 = t4.*u2.*u9;
t41 = -t49+u10;
t42 = t3.*t41.*u2;
t50 = t4.*u3.*u9;
t43 = -t50+u11;
t44 = t3.*t43.*u3;
t45 = t42+t44;
t46 = t5.*t11.*u9.*(1.0./2.0);
t52 = t10.*u9;
t53 = t45.*u1;
t47 = t46-t52-t53+u12;
t55 = t5.*t11.*(1.0./2.0);
t56 = t38.*u1;
t48 = t7+t9-t55-t56;
t51 = t2.*t14.*u9;
t54 = t51-t2.*t47.*u1;
fh = param6.*(u4-uh4)-t15.*t16.*t17.*t18.*(nl1.*param1.*t3.*t12.*t31+nl2.*param1.*t3.*t12.*t40.*t54);
if nargout > 1
    t57 = t2.*t14;
    t58 = t2.*t48.*u1;
    t59 = t57+t58;
    fh_udg = [t15.*t16.*t17.*t18.*(nl1.*param1.*t12.*t31.*t35.*2.0-nl1.*param1.*t3.*t12.*(t2.*t27+t2.*u1.*(t21+t23+u1.*(t6.*t39.*u5+t8.*t39.*u5-t20.*t35.*u2.*2.0-t22.*t35.*u3.*2.0)-t38.*u5+t5.*x1)-t2.*t48.*u5)+nl2.*param1.*t12.*t35.*t40.*t54.*2.0+nl2.*param1.*t3.*t12.*t40.*(t2.*t47+t2.*t48.*u9-t2.*u1.*(t42+t44+u1.*(t6.*t39.*u9+t8.*t39.*u9-t35.*t41.*u2.*2.0-t35.*t43.*u3.*2.0)-t38.*u9)));-t15.*t16.*t17.*t18.*(nl1.*param1.*t3.*t12.*(t2.*u1.*(u1.*(t3.*t20-t35.*u2.*u5)+t3.*u2.*u5)-t2.*t4.*u2.*u5)+nl2.*param1.*t3.*t12.*t40.*(t2.*u1.*(u1.*(t3.*t41-t35.*u2.*u9)+t3.*u2.*u9)-t2.*t4.*u2.*u9));-t15.*t16.*t17.*t18.*(nl1.*param1.*t3.*t12.*(t2.*u1.*(u1.*(t3.*t22-t35.*u3.*u5)+t3.*u3.*u5)-t2.*t4.*u3.*u5)+nl2.*param1.*t3.*t12.*t40.*(t2.*u1.*(u1.*(t3.*t43-t35.*u3.*u9)+t3.*u3.*u9)-t2.*t4.*u3.*u9));param6-t15.*t16.*t17.*t18.*(nl1.*param1.*t2.*t3.*t12.*u5+nl2.*param1.*t2.*t3.*t12.*t40.*u9);-nl1.*param1.*t3.*t15.*t16.*t18.*t59;-nl1.*param1.*t3.*t15.*t16.*u2;-nl1.*param1.*t3.*t15.*t16.*u3;nl1.*param1.*t4.*t15.*t16;-nl2.*param1.*t3.*t15.*t16.*t18.*t40.*t59;-nl2.*param1.*t3.*t15.*t16.*t40.*u2;-nl2.*param1.*t3.*t15.*t16.*t40.*u3;nl2.*param1.*t4.*t15.*t16.*t40];
end
if nargout > 2
    fh_uh = [zero;zero;zero;-one.*param6];
end
fh = reshape(fh,ng,1);
fh_udg = reshape(fh_udg,ng,1,nc);
fh_uh = reshape(fh_uh,ng,1,nch);
