function inelement(px, p)

S = [0 0; 1 0; 0 1; 1 1];
p1 = mapp(p,S);

% p2 = mapp(p,R);
% 
% 
% S = [0 0; 1 0; 0 1; 1 1];
% L = [-1 0; 0 0; -1 1; 0.5 1];
% R = [0 0; 1.0 0; 0.5 1; 1 1];
% xrect = [-1 1 0 1];
% 
% m = 11; n = 11; 
% [p,t] = squaremesh(m,n,0,1);
% % p(:,1) = xrect(1) + (xrect(2)-xrect(1))*p(:,1);
% % p(:,2) = xrect(3) + (xrect(4)-xrect(3))*p(:,2);
% 
% p1 = mapp(p,L);
% p2 = mapp(p,R);
