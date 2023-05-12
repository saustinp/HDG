
global Nc;
Nc = 80;    % Cut-off

global sigma_N;
sigma_N = randn(2*Nc+1,1);
sigma_N(Nc+1) = 0;

global deltaT;
deltaT = dt;