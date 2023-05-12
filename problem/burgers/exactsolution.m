function [up, um] = exactsolution(x, t, alpha, beta, gamma)

% alpha*t^2*u^2 + (- beta*t - 2*alpha*t*x - 1)*u + alpha*x^2 + beta*x + gamma = 0
a = alpha*t^2;
b = (- beta*t - 2*alpha*t*x - 1);
c = alpha*x.^2 + beta*x + gamma;

up = (-b + sqrt(b.*b - 4*a*c))/(2*a);
um = (-b - sqrt(b.*b - 4*a*c))/(2*a);





