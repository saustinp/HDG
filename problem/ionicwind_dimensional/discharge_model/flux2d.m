function [out,jac_out] = flux2d(p,udg,param,time)
%FLUX2D
%    [OUT,JAC_OUT] = FLUX2D(P,UDG,PARAM,TIME)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    20-Jul-2023 19:34:48
[ng,nc] = size(udg);
nch = 4;
nd = 2;
one = ones(ng,1);
param1 = param{1};
param2 = param{2};
param3 = param{3};
param4 = param{4};
param15 = param{15};
param18 = param{18};
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
x1 = p(:,1);
x3 = p(:,3);
x4 = p(:,4);
zero = zeros(ng,1);
t2 = one.*x1;
t3 = u8+x3;
t4 = u12+x4;
t5 = log(1.0e+1);
t12 = param1.*u3.*x1;
t13 = param2.*u2.*x1;
t14 = 1.0./param18;
t10 = param3.*t2;
t11 = param4.*t2;
t15 = -t2;
t16 = t3.^2;
t17 = t4.^2;
t18 = -t13;
t19 = 1.0./t5;
t22 = t16+t17;
t23 = 1.0./t22;
t24 = sqrt(t22);
t25 = param15.*t14.*t24.*1.0e+21;
t26 = log(t25);
t27 = t19.*t26.*1.0e+3;
t35 = t19.*t26.*3.477702502838188e+23;
t36 = t19.*t26.*4.633837111420517e+24;
t37 = t19.*t26.*3.4808749e-1;
t28 = t27-1.9e+3;
t38 = -t37;
t40 = t35+1.413709274499572e+24;
t41 = t36-6.700654851204797e+24;
t29 = tanh(t28);
t39 = 1.0e+1.^t38;
t30 = t29.^2;
t31 = t29./2.0;
t42 = t14.*t39.*u1.*5.849682301203802e+24;
t32 = t30-1.0;
t33 = t31+1.0./2.0;
t34 = t31-1.0./2.0;
t43 = t3.*t19.*t23.*t34.*3.477702502838188e+23;
t44 = t3.*t19.*t23.*t33.*4.633837111420517e+24;
t45 = t4.*t19.*t23.*t34.*3.477702502838188e+23;
t46 = t4.*t19.*t23.*t33.*4.633837111420517e+24;
t49 = t34.*t40;
t50 = t33.*t41;
t52 = t3.*t19.*t23.*t32.*t40.*5.0e+2;
t53 = t4.*t19.*t23.*t32.*t40.*5.0e+2;
t56 = t3.*t19.*t23.*t32.*t41.*5.0e+2;
t57 = t4.*t19.*t23.*t32.*t41.*5.0e+2;
t47 = -t44;
t48 = -t46;
t51 = -t50;
t54 = -t52;
t55 = -t53;
t58 = t49+t51;
out = [-x1.*(t14.*t58.*u5.*1.0e+1+t3.*t14.*t39.*u1.*5.849682301203802e+24);x1.*(param3.*u6-param2.*t3.*u2);x1.*(param4.*u7+param1.*t3.*u3);-u8.*x1;-x1.*(t14.*t58.*u9.*1.0e+1+t4.*t14.*t39.*u1.*5.849682301203802e+24);x1.*(param3.*u10-param2.*t4.*u2);x1.*(param4.*u11+param1.*t4.*u3);-u12.*x1];
if nargout > 1
    t61 = t43+t47+t54+t56;
    t62 = t45+t48+t55+t57;
    t59 = t14.*t58.*x1.*1.0e+1;
    t60 = -t59;
    mt1 = [t3.*t14.*t39.*x1.*-5.849682301203802e+24;zero;zero;zero;t4.*t14.*t39.*x1.*-5.849682301203802e+24;zero;zero;zero;zero;-param2.*t3.*x1;zero;zero;zero;-param2.*t4.*x1;zero;zero;zero;zero;param1.*t3.*x1;zero;zero;zero;param1.*t4.*x1;zero;zero;zero;zero;zero;zero;zero;zero;zero;t60;zero;zero;zero;zero;zero;zero;zero;zero;t10;zero;zero;zero;zero;zero;zero;zero;zero;t11;zero;zero;zero;zero;zero;-x1.*(t42+t14.*t61.*u5.*1.0e+1-t14.*t16.*t23.*t39.*u1.*2.036201229523455e+24);t18;t12;t15];
    mt2 = [-x1.*(t14.*t61.*u9.*1.0e+1-t3.*t4.*t14.*t23.*t39.*u1.*2.036201229523455e+24);zero;zero;zero;zero;zero;zero;zero;t60;zero;zero;zero;zero;zero;zero;zero;zero;t10;zero;zero;zero;zero;zero;zero;zero;zero;t11;zero;-x1.*(t14.*t62.*u5.*1.0e+1-t3.*t4.*t14.*t23.*t39.*u1.*2.036201229523455e+24);zero;zero;zero];
    mt3 = [-x1.*(t42+t14.*t62.*u9.*1.0e+1-t14.*t17.*t23.*t39.*u1.*2.036201229523455e+24);t18;t12;t15];
    jac_out = [mt1;mt2;mt3];
end

out = reshape(out,ng,nch,nd);
jac_out = reshape(jac_out,ng,nch,nd,nc);
end
