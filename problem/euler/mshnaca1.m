function mesh=mshnaca1(angle,hs,porder,xd,thick)


if nargin<1, angle=0; end
if nargin<2, hs=[0.009 0.04 0.08 0.35 0.5 1]; end
if nargin<3, porder=1; end
if nargin<4, xd=[-4, 5, -4, 4]; end
if nargin<5, thick=12; end

xlim=[0,1.0089304129];
%xlim=[0,1];
a=(thick/100)/.2*[0.2969,-0.1260,-0.3516,0.2843,-0.1015];
theta=angle*pi/180;

x1=(hs(1)/a(1))^2;
x2=xlim(2)-hs(2);
x=[xlim(1),x1,x2,xlim(2)]';
[x,y]=limit_arcratio(x,hs(3),a);
R=[cos(theta),-sin(theta); sin(theta),cos(theta)];

pv1=[x(end:-1:1),y(end:-1:1); x(2:end-1),-y(2:end-1)]*R;
pv2=[xd(1),xd(3); xd(2),xd(3); xd(2),xd(4); xd(1),xd(4)];

[mesh.p,mesh.t]=polymesh({pv1,pv2},[1,1],[0,1;1,0],[0.85,1.4],@href,hs);

[mesh.p,mesh.t] = fixmesh(mesh.p,mesh.t);
[mesh.f,mesh.t2f] = mkt2f(mesh.t);

bndexpr = {'all(sqrt((p(:,1)-0.5).^2+p(:,2).^2)<0.6)', ...
           'all(sqrt((p(:,1)-0.5).^2+p(:,2).^2)>0.6)'};     
mesh.f = setbndnbrs(mesh.p,mesh.f,bndexpr);

mesh.fcurved = (mesh.f(:,4)<0);
ic = mesh.fcurved;
mesh.tcurved = repmat(false,size(mesh.t,1),1);
mesh.tcurved(mesh.f(ic,3)) = true;

mesh.porder = porder;
[mesh.plocal,mesh.tlocal] = uniformlocalpnts(mesh.porder);
mesh.dgnodes = createnodes(mesh,@fdnaca,xd,thick);
mesh.elemtype = 0;

% function d=dnaca(p,a)
% 
% e1=polyval([a(5:-1:2),0],p(:,1));
% e2=a(1)^2*p(:,1);
% d1=(p(:,2)-e1).^2-e2;
% d2=(p(:,2)+e1).^2-e2;
% d=max(d1,d2);


function d=fdnaca(p,xd,thick)
%FDNACA Distance function for NACA foil inside a rectangular domain. 
%  D=FDNACA(P,XD,THICK)
%
%      P:         Node positions (Nx2)
%      XD(4):     [x_min, x_max, y_min, y_max] for far field boundary
%      THICK:     Airfoil thickness in percentage
%
%   See also:  NACA, DPOLY, DRECTANGLE, DDIFF
%
th=pi:-pi/200:pi/2;
xt = cos(th)+1; xt(end)=[];  yt=naca(xt,thick);  
xb=fliplr(xt); yb=-naca(xb,thick);
pv=[xt',yt';1.0089304129,0;xb',yb'];
dfoil=dpoly(p,pv);                          % distance to foil

drec=drectangle(p,xd(1),xd(2),xd(3),xd(4));     % distance to domain edge

d=ddiff(drec,dfoil);

%--------------------------------------------------------------------------

function y=naca(x,thick)
%NACA Equation for a NACA profile. 
%  Y=NACA(X,THICK)
%
%      Y:         y-coordinate
%      X:         x-coordinate
%      THICK:     Airfoil thickness in percentage
%
y=5*0.01*thick*(0.29690*sqrt(x)-0.12600*x-0.35160*x.^2+0.28430*x.^3-0.10150*x.^4);

%--------------------------------------------------------------------------

function h=href(p,hs)



pv1=[3,.75;0,.25;0,-.25;3,-.75;3,.75];

pv2=[2,.5;0,.25;0,-.25;2,-.5;2,.5];

inside1=inpolygon(p(:,1),p(:,2),pv1(:,1),pv1(:,2));

inside2=inpolygon(p(:,1),p(:,2),pv2(:,1),pv2(:,2));

h=inf*ones(size(p,1),1);

h(inside1)=hs(5);

h(inside2)=hs(4);



function [x,y]=limit_arcratio(x,hmax,a)



while 1

  y=polyval([a(5:-1:2),0],x)+a(1)*sqrt(x);
  %y=naca(x,12);  
  s=sqrt((x(2:end)-x(1:end-1)).^2+(y(2:end)-y(1:end-1)).^2);

  

  [maxs,ix]=max(s);

  if maxs>1.4*hmax

    x=[x(1:ix);(x(ix)+x(ix+1))/2;x(ix+1:end)];

    continue;

  end

  

  ratio=s(2:end)./s(1:end-1);

  ix=find(ratio>2 | ratio<1/2);

  

  if isempty(ix)

    break;

  else

    ix=ix(1);

    if ratio(ix)>2

      x=[x(1:ix+1);(x(ix+1)+x(ix+2))/2;x(ix+2:end)];

    else

      x=[x(1:ix);(x(ix)+x(ix+1))/2;x(ix+1:end)];

    end

  end

end

