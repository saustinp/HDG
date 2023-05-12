function [] = test_mFlux()

g = 1;
e = 2;
c = 3;
d = 4;

n(g) = 3;
n(e) = 3;
n(c) = 2;
n(d) = 2;

A = mArray.ones(n,[g e c]);

B = mFunction(@foo,n,A,c,[c d]);

B

function [B] = mFunction(f,n,A,cc,dd)

% B(dd(1),1,dd(2),1) = 1.*A(cc,1);
% B(dd(1),1,dd(2),2) = 2.*A(cc,1);
% 
% B(dd(1),2,dd(2),1) = 3.*A(cc,2);
% B(dd(1),2,dd(2),2) = 4.*A(cc,2);

B = feval(f,n,A,cc,dd);

function b = foo(n,a,cc,dd)

B = mArray.zerosTransform(n,cc,dd);

c = cc(1);
d = dd(2);

rho  = a(c,1);
rhou = a(c,2);

x = 1;
y = 2;

f1x = 1.*rho;
f1y = 2.*rho;
f2x = 3.*rhou;
f2y = 4.*rhou;

b(c,rho,d,x) = f1x;
b(c,rho,d,y) = f1y;
b(c,rhou,d,x) = f2x;
b(c,hrou,d,y) = f2y;

