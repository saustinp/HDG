load('/Users/peraire/Dropbox (MIT)/DATA/Buffeting Foil/foilgeom.mat')
clf;

xf =[xf, yf];
clear yf;
t = 0:length(xf(:,1))-1;
spf = [spline(t,xf(:,1)),spline(t,xf(:,2))];
spfder = [fnder(spf(1)),fnder(spf(2))];
dxf = [fnval(spfder(1),t)', fnval(spfder(2),t)'];
dsf = sqrt(dxf(:,1).^2+dxf(:,2).^2);
dxf = [dxf(:,1)./dsf, dxf(:,2)./dsf];
nf = -[dxf(:,2), -dxf(:,1)];

plot(xf(:,1),xf(:,2),'-o',xf(:,1)+0.1*nf(:,1),xf(:,2)+0.1*nf(:,2),'-o');
hold on;

%BOTTOM WAKE
tw = 0:2;
xwb = zeros(3,2);
xwb(1,:) = xf(1,:);
xwb(2,:) = [xf(1,1) + 0.3, xf(1,2)-0.01];
xwb(3,:) = [xf(1,1) + 1.0, xf(1,2)+0.0];
spwb = [spline(tw,xwb(:,1)),spline(tw,xwb(:,2))];
spwbder = [fnder(spwb(1)),fnder(spwb(2))];

twp = 0:0.1:2;
xwbp = [fnval(spwb(1),twp)', fnval(spwb(2),twp)'];
dxwb = [fnval(spwbder(1),twp)', fnval(spwbder(2),twp)'];
dswb = sqrt(dxwb(:,1).^2+dxwb(:,2).^2);
nb = [dxwb(:,2)./dswb, -dxwb(:,1)./dswb];
nb = ([nb(:,1).*twp', nb(:,2).*twp']/2 + [nf(1,1)*(2-twp'), nf(1,2)*(2-twp')]);
nbs = sqrt(nb(:,1).^2+nb(:,2).^2);
nb = 0.1*nb./nbs;
plot(xwbp(:,1),xwbp(:,2),'-+',xwbp(:,1)+nb(:,1),xwbp(:,2)+nb(:,2),'-+')

for i = 1:length(twp)
     vectarrow([xwbp(i,1),xwbp(i,2)],[xwbp(i,1)+nb(i,1),xwbp(i,2)+nb(i,2)]);
     hold on;
end

%TOP WAKE
xwt = zeros(3,2);
xwt(1,:) = xf(end,:);
xwt(2,:) = [xf(end,1) + 0.3, xf(end,2)-0.0];
xwt(3,:) = [xf(end,1) + 1.0, xf(end,2)+0.03];
spwt = [spline(tw,xwt(:,1)),spline(tw,xwt(:,2))];
spwtder = [fnder(spwt(1)),fnder(spwt(2))];

xwtp = [fnval(spwt(1),twp)', fnval(spwt(2),twp)'];
dxwt = [fnval(spwtder(1),twp)', fnval(spwtder(2),twp)'];
dswt = sqrt(dxwt(:,1).^2+dxwt(:,2).^2);
nt = -[dxwt(:,2)./dswt, -dxwt(:,1)./dswt];
nt = ([nt(:,1).*twp', nt(:,2).*twp']/2 + [nf(end,1)*(2-twp'), nf(end,2)*(2-twp')]);
nts = sqrt(nt(:,1).^2+nt(:,2).^2);
nt = 0.1*nt./nts;

plot(xwtp(:,1),xwtp(:,2),'-+',xwtp(:,1)+nt(:,1),xwtp(:,2)+nt(:,2),'-+')

for i = 1:length(twp)
     vectarrow([xwtp(i,1),xwtp(i,2)],[xwtp(i,1)+nt(i,1),xwtp(i,2)+nt(i,2)]);
     hold on;
end

nw = 20;
ttwp = 0:20;
d0 = 0.025;
rat = 1.05;
ttwp = d0*(rat.^ttwp -1)/(rat-1);
ttwp = 2*ttwp/ttwp(end);

xwttp = [fnval(spwt(1),ttwp)', fnval(spwt(2),ttwp)'];
plot(xwttp(:,1),xwttp(:,2),'*');
xwtbp = [fnval(spwb(1),ttwp)', fnval(spwb(2),ttwp)'];
plot(xwtbp(:,1),xwtbp(:,2),'*');


axis equal;
grid on;

% xwt = xf(1,:) - xw(1,:);
% xwb = xf(end,:) - xw(1,:);
% 
% spwx = spline(tw,xw(:,1));
% spwy = spline(tw,xw(:,2));
% spwxder = fnder(spwx);
% spwyder = fnder(spwy);
% 
% xwp = zeros(length(ttwp),2);
% xwp(:,1) = fnval(spwx,ttwp);
% xwp(:,2) = fnval(spwy,ttwp);
% dwx = fnval(spwxder,ttwp);
% dwy = fnval(spwyder,ttwp);
% dws = sqrt(dwx.^2 + dwy.^2);
% dwx = dwx./dws;
% dwy = dwy./dws;
% xwt = xwp + xwt;
% xwb = xwp + xwb;
% xwot = xwt + 0.1*[-dwy' + dwx'];
% xwob = xwt - 0.1*[-dwy' + dwx'];
% plot(xwot(:,1),xwot(:,2),'-o',xwob(:,1),xwob(:,2),'-o');
% 
% 
% 
% 
% ttp = distribute(80,spx,spy,n);
% xp = fnval(spx,ttp);
% yp = fnval(spy,ttp);
% 
% plot(xp,yp,'*');
% hold on;
% axis equal;
% 
% dx = fnval(spxder,t);
% dy = fnval(spyder,t);
% ds = sqrt(dx.^2 + dy.^2);
% dx = dx./ds;
% dy = dy./ds;
% xfo = xf(:,1)-0.1*dy';
% yfo = xf(:,2)+0.1*dx';
% 
% spxo = spline(t,xfo);
% spyo = spline(t,yfo);
% xpo = fnval(spxo,t);
% ypo = fnval(spyo,t);
% plot(xpo,ypo,'o');
% 
% for i = 1:length(t)
%     vectarrow([xf(i),yf(i)],[xfo(i),yfo(i)]);
%     hold on;
% end
% 
% axis equal;