function f = limiting2(x,x1,x2,fmin,fmax)

m = (fmax-fmin)/(x2-x1);
alpha = 100;

f = m.*(x-x1).*(atan(alpha.*(x-x1))/pi +.5) + m*(-atan(alpha)/pi +.5) + fmin;
f = lmin(f-fmax,100*alpha) + fmax;

% To test:
% x=.1e6:.1e6:5e7;
% f = limiting2(x, 2e6, 2e7,0, .05);
% plot(x,f)