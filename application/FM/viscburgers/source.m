function [sr,sr_udg] = source(p,udg,param,time)

[ng,nc] = size(udg);
nch = 1;

kappa = param{1};

x = p(:,1);
y = p(:,2);

%sr = x.*y.^2.*tanh((x - 1)/kappa).^2.*tanh((y - 1)/kappa).^2 - x.*tanh((x - 1)/kappa).*((2.*y.*tanh((y - 1)/kappa).*(tanh((y - 1)/kappa).^2 - 1))/kappa - 2.*tanh((y - 1)/kappa).^2 + 2) - y.*tanh((y - 1)/kappa).*((2.*x.*tanh((x - 1)/kappa).*(tanh((x - 1)/kappa).^2 - 1))/kappa - 2.*tanh((x - 1)/kappa).^2 + 2) + x.^2.*y.*tanh((x - 1)/kappa).^2.*tanh((y - 1)/kappa).^2 - (x.^2.*y.^2.*tanh((x - 1)/kappa).*tanh((y - 1)/kappa).^2.*(tanh((x - 1)/kappa).^2 - 1))/kappa - (x.^2.*y.^2.*tanh((x - 1)/kappa).^2.*tanh((y - 1)/kappa).*(tanh((y - 1)/kappa).^2 - 1))/kappa;

sr = ones(ng,nch);

sr_udg = zeros(ng,nch,nc); 




