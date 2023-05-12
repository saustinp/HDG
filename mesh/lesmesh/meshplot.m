function hh=meshplot(p,t)
%MESHPLOT  Plot Mesh Structure (with straight edges)
%    MESHPLOT(P,T)

hh=[];
pars={'facecolor',[.8,1,.8],'edgecolor','k','Linew',1};
clf,hh=[hh;patch('faces',t,'vertices',p,pars{:})];
view(2),axis equal

% if opts(1)
%   xx=squeeze(mesh.dgnodes(:,1,:));
%   yy=squeeze(mesh.dgnodes(:,2,:));
%   zz=0*xx;
%   line(xx(:),yy(:),zz(:),'lines','n','marker','.','markersi',16,'col','b');
% end

