function [p,t] = lesmesh
mw = 30;
mf = 80;
m = mf+2*mw;
n = 20;

load('foilgeom.mat')

ns = length(xf);
t = 0:ns-1;
spx = spline(t,xf);
spy = spline(t,yf);

% plot(xf,yf);
% hold on;
% axis equal;

[Y,X] = meshgrid(0:1/n:1,0:1/m:1);

%
% Uniform distribution over foil
%
ttp = distribute(mf,spx,spy,ns);
xv = zeros(mf+2*mw+1,1);
xv(mw+1:mw+mf+1,1) = 2*ttp/(ns-1);

%
% Geometric distribution on wakes
%

ttwp = 0:mw;
d0 = 0.025;
rat = 1.05;
ttwp = d0*(rat.^ttwp -1)/(rat-1);
ttwp = ttwp/ttwp(end);
xv(1:mw+1,1) = -fliplr(ttwp);
xv(mw+mf+1:end,1) = ttwp + 2;

X = xv*ones(1,n+1);

yv = 0:n;
dlay = 0.1;
rat = 1.3;
d0 = 1e-4;
yvu = yv*dlay/yv(end);
yv = d0*(rat.^yv-1)/(rat-1);
yv = yv*dlay/yv(end);

wg = ones(mf+2*mw+1,1);
wg(1:mw+1) = 0:1/mw:1;
wg(mw+mf+1:end) = 1:-1/mw:0;
wg1 = ones(mf+2*mw+1,1)-wg;
al = 0.3;

Y = (wg+al*wg1)*yv + (1-al)*wg1*yvu;

p = [X(:),Y(:)];


dlay = 0.1;
rat = 1.3;
d0 = 1e-4;
p(:,2) = p(:,2)*(n+1);
p(:,2) = d0*(rat.^p(:,2)-1)/(rat-1);
p(:,2) = p(:,2)*dlay/p(end,2);

% pp(:,2) = p(:,2)*dlay/p(end,2);

t = [1, 2, m+3, m+2];
t = kron(t,ones(m,1)) + kron(ones(size(t)),(0:m-1)');
t = kron(t,ones(n,1)) + kron(ones(size(t)),(0:n-1)'*(m+1));
level = kron((n:-1:1)',ones(1,m)); level = level(:);

t2t = mkt2t(t,1);
% meshplot(p,t);
% grid on;



%
% Trailing edge refinement
%
for i = 1:2
    pe = 0.25*[p(t(:,1),1) + p(t(:,2),1) + p(t(:,3),1) + p(t(:,4),1),  ...
               p(t(:,1),2) + p(t(:,2),2) + p(t(:,3),2) + p(t(:,4),2)];
    ind = 1:size(t,1);
    il = pe(:,1) > 0;
    indl = ind(il);
    [ymx,ie] = max(pe(il,2)-pe(il,1));
    [p,t,t2t,level] = refineall(indl(ie),p,t,t2t,level);
    
    pe = 0.25*[p(t(:,1),1) + p(t(:,2),1) + p(t(:,3),1) + p(t(:,4),1),  ...
               p(t(:,1),2) + p(t(:,2),2) + p(t(:,3),2) + p(t(:,4),2)];
    ind = 1:size(t,1);
    il = pe(:,1) < 2;
    indl = ind(il);
    [ypx,ie] = max(pe(il,2)+pe(il,1));
    [p,t,t2t,level] = refineall(indl(ie),p,t,t2t,level);
end

%
% Shock Refinement
%

nt = size(t,1);
pe = 0.25*[p(t(:,1),1) + p(t(:,2),1) + p(t(:,3),1) + p(t(:,4),1),  ...
           p(t(:,1),2) + p(t(:,2),2) + p(t(:,3),2) + p(t(:,4),2)];
for ie = 1:nt
    if (t2t(ie,3) == 0),
        if pe(ie,1) > 1.48 & pe(ie,1) < 2,
            [p,t,t2t,level] = refineall(ie,p,t,t2t,level);
        end
    end
end 

nt = size(t,1);
pe = 0.25*[p(t(:,1),1) + p(t(:,2),1) + p(t(:,3),1) + p(t(:,4),1),  ...
           p(t(:,1),2) + p(t(:,2),2) + p(t(:,3),2) + p(t(:,4),2)];
for ie = 1:nt
    if (t2t(ie,3) == 0),
        if pe(ie,1) > 1.5 & pe(ie,1) < 1.8,
            [p,t,t2t,level] = refineall(ie,p,t,t2t,level);
        end
    end
end  

nt = size(t,1);
for ie = 1:nt
    if (level(ie) == 4),
            [p,t,t2t,level] = refine(ie,p,t,t2t,level);
    end
end  


pe = 0.25*[p(t(:,1),1) + p(t(:,2),1) + p(t(:,3),1) + p(t(:,4),1),  ...
           p(t(:,1),2) + p(t(:,2),2) + p(t(:,3),2) + p(t(:,4),2)];     
ind = 1:size(t,1);
il = (level == 8 & pe(:,1) > 2.0);
indl = ind(il);
[ymx,ie] = min(pe(il,1));
[p,t,t2t,level] = refine(indl(ie),p,t,t2t,level);

pe = 0.25*[p(t(:,1),1) + p(t(:,2),1) + p(t(:,3),1) + p(t(:,4),1),  ...
           p(t(:,1),2) + p(t(:,2),2) + p(t(:,3),2) + p(t(:,4),2)];
ind = 1:size(t,1);
il = (level == 8 & pe(:,1) < 0.0);
indl = ind(il);
[ymx,ie] = max(pe(il,1));
[p,t,t2t,level] = refine(indl(ie),p,t,t2t,level);

pe = 0.25*[p(t(:,1),1) + p(t(:,2),1) + p(t(:,3),1) + p(t(:,4),1),  ...
           p(t(:,1),2) + p(t(:,2),2) + p(t(:,3),2) + p(t(:,4),2)];
ind = 1:size(t,1);
il = (level == 8 & pe(:,1) > 0.0);
indl = ind(il);
[ymx,ie] = min(pe(il,1));
[p,t,t2t,level] = refine(indl(ie),p,t,t2t,level);

pe = 0.25*[p(t(:,1),1) + p(t(:,2),1) + p(t(:,3),1) + p(t(:,4),1),  ...
           p(t(:,1),2) + p(t(:,2),2) + p(t(:,3),2) + p(t(:,4),2)];
nt = size(t,1);       
for ie = 1:nt
    if (level(ie) == 8),
        if pe(ie,1) > 0.7 & pe(ie,1) < 2.0,
            [p,t,t2t,level] = refine(ie,p,t,t2t,level);
        end
    end
end  
% 
% meshplot(p,t);
% grid on;
% pause;


%
% Now the trailing edge mesh
%

ilb = (p(:,2) == 0 & p(:,1) <= 0);
ilt = (p(:,2) == 0 & p(:,1) >= 2);

p = map(p);

% meshplot(p,t);
% grid on;

pt = p(ilt,:);
[dum,ind] = sort(pt,1);
pt = pt(ind(:,1),:);
pb = p(ilb,:);
[dum,ind] = sort(pb,1);
pb = pb(ind(:,1),:);
nx = size(pt,1)-1;

nlw = 7;
pw = pb;
for i = 1:nlw
    pw = [pw; i*pt/nlw + (nlw-i)*pb/nlw];
end

tw = [1, 2, nx+3, nx+2];
tw = kron(tw,ones(nx,1)) + kron(ones(size(tw)),(0:nx-1)');
tw = kron(tw,ones(nlw,1)) + kron(ones(size(tw)),(0:nlw-1)'*(nx+1));

np = size(p,1);
p = [p; pw];
t = [t; tw+np];

% meshplot(p,t);
% grid on;
% size(p)
% size(t)

[p,t] = fixmesh(p,t);
% meshplot(p,t);
% grid on;
% size(p)
% size(t)

% [p,t,t2t,level,pp] = refine(5,p,t,t2t,level,pp);
% [p,t,t2t,level,pp] = refine(9,p,t,t2t,level,pp);
% [p,t,t2t,level,pp] = refine(13,p,t,t2t,level,pp);
% [p,t,t2t,level,pp] = refine(17,p,t,t2t,level,pp);
% 
% % meshplot(p,t);
% % axis equal
% 
% %
% % Map foil
% %
% spxder = fnder(spx);
% spyder = fnder(spy);
% 
% pm = 0*p;
% pm(:,1) = fnval(spx,pp(:,1));
% pm(:,2) = fnval(spy,pp(:,1));
% dx = fnval(spxder,pp(:,1));
% dy = fnval(spyder,pp(:,1));
% ds = sqrt(dx.^2 + dy.^2);
% dx = dx./ds;
% dy = dy./ds;
% pm(:,1) = pm(:,1)-pp(:,2).*dy;
% pm(:,2) = pm(:,2)+pp(:,2).*dx;
% 
% meshplot(pm,t);
% axis equal