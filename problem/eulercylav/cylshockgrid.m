function [pn,tn,dgnodes] = cylshockgrid(tm,ds,polyn,porder)

nx = cos(tm);
ny = sin(tm);
dm = polyval(polyn, tm); 
xm = dm.*cos(tm);
ym = dm.*sin(tm);
xn = zeros(length(ds), length(tm));
yn = 0*xn;
for i = 1:length(ds)
    xn(i,:) = xm + ds(i)*nx;
    yn(i,:) = ym + ds(i)*ny;
end
[~,tn] = squaremesh(length(ds),length(tm),0,1);
pn = [xn(:) yn(:)];

[~, t0] =masternodes(length(tm)-1,1);
mesh0 = mkmesh(tm,t0,porder,{'true'},1,1);
dgt = mesh0.dgnodes(:);

[~, t0] =masternodes(length(ds)-1,1);
mesh0 = mkmesh(ds(:),t0,porder,{'true'},1,1);
dgs = mesh0.dgnodes(:);

nx = cos(dgt);
ny = sin(dgt);
dm = polyval(polyn, dgt); 
xm = dm.*cos(dgt);
ym = dm.*sin(dgt);
xn = zeros(length(dgs), length(dgt));
yn = 0*xn;
for i = 1:length(dgs)
    xn(i,:) = xm + dgs(i)*nx;
    yn(i,:) = ym + dgs(i)*ny;
end
p1 = porder+1;
x = reshape(xn,[p1 length(ds)-1 p1 length(tm)-1]);
x = reshape(permute(x,[1 3 2 4]),[p1*p1 (length(tm)-1)*(length(ds)-1)]);
y = reshape(yn,[p1 length(ds)-1 p1 length(tm)-1]);
y = reshape(permute(y,[1 3 2 4]),[p1*p1 (length(tm)-1)*(length(ds)-1)]);
dgnodes = zeros(p1*p1,2,(length(tm)-1)*(length(ds)-1));
dgnodes(:,1,:) = x;
dgnodes(:,2,:) = y;

% figure(2); clf; hold on;
% simpplot(pn,tn);
% plot(x(:),y(:),'o');
% axis equal; axis tight;

%mesh = mkmesh(pn,tn,porder,{'true'},1,1);
%[mesh.dgnodes(:,:,10)-dgnodes(:,:,10)]
