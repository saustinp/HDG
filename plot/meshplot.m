function hh = meshplot(mesh,opts,udg,elemsToPlot)
%MESHPLOT  Plot mesh structure 
%    MESHPLOT(MESH,[OPTS])
%
%    MESH:       Mesh structure
%    OPTS:       (logical)
%      OPTS(1):  Plot elements/faces using p (default) or dgnodes
%      OPTS(2):  Plot dgnodes 
%      OPTS(3):  Plot element numbers (2D only)
%      OPTS(4):  Plot node numbers (2D only)
%      OPTS(5):  Plot face numbers (2D only)
% elemsToPlot: vector of element indices to highlight in red

if nargin<2 || isempty(opts), opts=0; end
if length(opts)<2, opts=[opts,0]; end
if length(opts)<3, opts=[opts,0]; end
if length(opts)<4, opts=[opts,0]; end
if length(opts)<5, opts=[opts,0]; end

p=mesh.p;
t=mesh.t;
dim=size(p,2);
dpl=size(mesh.plocal,2);

% surface mesh
if dpl==2 && dim==3
    surface = 1;
    f=mesh.f;
else 
    surface = 0;
end
    
if dim < 1 || dim > 3,
    error('Only can handle dim=1, dim=2 or dim=3');
end

%pars={'facecolor',[0.2,0.5,1.0],'edgecolor','r','Linew',0.5};
pars={'facecolor',[1,1,1],'edgecolor','k','Linew',0.5};

if exist('hh','var')==0
    hh=[];
end
if dim==1
    plot(p,0*p);
elseif dim == 2 || surface==1    
    if opts(1) == 0 % plot elements using p        
        hh=[hh;patch('faces',t,'vertices',p,pars{:})];
        pars={'facecolor','r','edgecolor','r','Linew',0.5};
        t = t(elemsToPlot,:);
        hh=[[];patch('faces',t,'vertices',p,pars{:})];    
        
    else            % plot elements using dgnodes
        e=boundedges(mesh.plocal,mesh.tlocal,mesh.elemtype);
        e1=segcollect(e);        
        axis equal,axis off
        nt=size(mesh.dgnodes,3);
        hh=zeros(nt,1);
        for it=1:nt
            px=mesh.dgnodes(:,1,it);
            py=mesh.dgnodes(:,2,it);
            if surface==1  
                pz=mesh.dgnodes(:,3,it);
            else
                pz=0*px;
            end
            hh(it)=patch(px(e1{1}'),py(e1{1}'),pz(e1{1}'),0.0*e1{1}',pars{:});
        end        
    end
    if surface==0
        view(2),axis equal;
    else
        view(3),axis equal;
    end    
elseif dim ==3    
    if opts(1)==0 % plot boundary faces using p
        %f = mesh.f;
        clf,hh=[hh;patch('faces',mesh.f(:,1:end-2),'vertices',p,pars{:})];                            
    else          % plot boundary faces using dgnodes
        bf = find(mesh.f(:,end)<0);    
        e=boundedges(mesh.plocfc,mesh.tlocfc,mesh.elemtype);        
        e1=segcollect(e);        
        axis equal,axis off        
        nt=length(bf);
        hh=zeros(nt,1);
        for it=1:nt
            el = mesh.f(bf(it),end-1);
            fc = mesh.t2f(el,:);
            fi = find(fc==bf(it));
            px=mesh.dgnodes(mesh.perm(:,fi),1,el);
            py=mesh.dgnodes(mesh.perm(:,fi),2,el);
            pz=mesh.dgnodes(mesh.perm(:,fi),3,el);
            %pw=udg(mesh.perm(:,fi),1,el);
            %hh(it)=patch(px(e1{1}'),py(e1{1}'),pz(e1{1}'),pw(e1{1}'),pars{:});
            %hh(it)=patch(px,py,pz,0.0*px,pars{:});
            hh(it)=patch(px(e1{1}'),py(e1{1}'),pz(e1{1}'),0*pz(e1{1}'),pars{:});
        end
    end        
    view(3),axis equal; 
end

if opts(2)
    if dim ==1
        xx=squeeze(mesh.dgnodes(:,1,:));
        yy=0*xx;
        zz=0*xx;
    elseif dim == 2       
        xx=squeeze(mesh.dgnodes(:,1,:));
        yy=squeeze(mesh.dgnodes(:,2,:));
        zz=0*xx;
    elseif dim ==3
        xx=squeeze(mesh.dgnodes(:,1,:));
        yy=squeeze(mesh.dgnodes(:,2,:));
        zz=squeeze(mesh.dgnodes(:,3,:));
    end    
    line(xx(:),yy(:),zz(:),'lines','n','marker','.','markersize',16,'col','b');
end

if opts(3)
    if dim == 2 || surface==1   
        for it=1:size(t,1)
            pmid=mean(p(t(it,:),:),1);
            txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
            text(pmid(1),pmid(2),num2str(it),txtpars{:});
        end
    end
end

if opts(4)
    if dim == 2 || surface==1
        for i=1:size(p,1)
            txtpars={'fontname','times','fontsize',20,'fontweight','bold', ...
                'horizontala','center','col','w', 'BackgroundColor',[0.5,0.5,0.5]};
            text(p(i,1),p(i,2),num2str(i),txtpars{:});
        end
    end
end

if opts(5)
    if dim == 2 || surface==1
        for i=1:size(mesh.f,1)
            txtpars={'fontname','times','fontsize',20,'fontweight','bold', ...
                'horizontala','center','col','w', 'BackgroundColor',[0.5,0.5,0.5]};            
            fi = mesh.f(i,1:2);
            pm = (mesh.p(fi(1),:)+mesh.p(fi(2),:))/2;
            text(pm(1),pm(2),num2str(i),txtpars{:});
        end
    end
end

if nargout<1, clear hh; end
set(gcf,'color','w');


function e=boundedges(p,t,elemtype)
%BOUNDEDGES Find boundary edges from triangular mesh
%   E=BOUNDEDGES(P,T)

% Form all edges, non-duplicates are boundary edges

if elemtype==0
    edges=[t(:,[1,2]);
           t(:,[1,3]);
           t(:,[2,3])];
    node3=[t(:,3);t(:,2);t(:,1)];
else
    edges=[t(:,[1,2]);
           t(:,[2,3]);
           t(:,[3,4]);
           t(:,[4,1]);];
    node3=[t(:,4);t(:,3);t(:,2);t(:,1)];    
end
edges=sort(edges,2);
[foo,ix,jx]=unique(edges,'rows');
vec=histc(jx,1:max(jx));
qx=find(vec==1);
e=edges(ix(qx),:);
node3=node3(ix(qx));

% Orientation
v1=p(e(:,2),:)-p(e(:,1),:);
v2=p(node3,:)-p(e(:,1),:);
ix=find(v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1)>0);
e(ix,[1,2])=e(ix,[2,1]);


function e1=segcollect(e)
%SEGCOLLECT Collect polygons from edge segments.

ue=unique(e(:));
he=histc(e(:),ue);
current=ue(min(find(he==1))); % Find an endpoint
if isempty(current) % Closed curve
  current=e(1,1);
end
e1=current;
while ~isempty(e)
  ix=min(find(e(:,1)==e1(end)));
  if isempty(ix)
    ix=min(find(e(:,2)==e1(end)));
    if isempty(ix) % >1 disjoint curves, recur
      rest=segcollect(e);
      e1={e1,rest{:}};
      return;
    end
    next=e(ix,1);
  else
    next=e(ix,2);
  end
  e1=[e1,next];
  e(ix,:)=[];
end
e1={e1};




