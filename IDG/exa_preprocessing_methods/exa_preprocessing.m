function [app,mesh,master,dmd] = exa_preprocessing(app,mesh,master)
% master=[];
app.nd  = mesh.nd;
app.ncx = app.nd;
% [app.ne, app.nve] = size(mesh.t);
app.ne = mesh.ne;
app.nve = mesh.nve;

app.elemtype = mesh.elemtype;
% if (app.nd==2) && (app.nve==4)
%     app.elemtype=1;    
% end
% if (app.nd==3) && (app.nve==8)
%     app.elemtype=1;    
% end
app.porder = mesh.porder;
app.pgauss = 2*app.porder;

% % master struct
% master = Master(app);
% writemaster(master,filemaster,'native');        

% obtain the PDE model
% pdemodel = str2func(app.modelfile);
% pde = pdemodel();

app.boundaryconditions = app.bcm;
app.uinf = app.ui;
nuinf = length(app.uinf);
% nparam = length(app.physicsparam);
nparam = length(app.arg);
% xdgsym = sym('xdg',[app.ncx 1]); 
% uinfsym = sym('uinf',[nuinf 1]); 
% paramsym = sym('param',[nparam 1]); 
% if isfield(pde, 'initu')    
%     udgsym = pde.initu(xdgsym, paramsym, uinfsym);         
%     app.ncu = length(udgsym(:));
% else
%     error("pde.initu is not defined");
% end
% if isfield(pde, 'initv')
%     odgsym = pde.initv(xdgsym, paramsym, uinfsym);         
%     app.nco = length(odgsym(:));
% elseif isfield(mesh, 'vdg')            
%     app.nco = size(mesh.vdg,2);    
% else    
%     app.nco = 0;
% end
% if isfield(pde, 'initw')
%     wdgsym = pde.initw(xdgsym, paramsym, uinfsym);         
%     app.ncw = length(wdgsym(:));
% elseif isfield(mesh, 'wdg')            
%     app.ncw = size(mesh.wdg,2);        
% else    
%     app.ncw = 0;
% end

% if app.model=="ModelC" || app.model=="modelC"
%     app.wave = 0;
%     app.nc = app.ncu;
% elseif app.model=="ModelD" || app.model=="modelD"     
%     app.wave = 0;
%     app.nc = (app.ncu)*(app.nd+1);
% elseif app.model=="ModelW" || app.model=="modelW"
%     app.tdep = 1;
%     app.wave = 1;
%     app.nc = (app.ncu)*(app.nd+1);
% end
app.ncq = app.nc - app.ncu;
app.nch  = app.ncu;               

% if max(app.dt)>0
%     app.tdep = 1;
% else
%     app.tdep = 0;
% end

disp('run facenumbering...');  
disp("TODO: no periodic BCs")
[ftmp, mesh.tprd, t2t] = facenumbering(mesh.p',mesh.t',app.elemtype,mesh.bndexpr,{});

mpiprocs = 1;
%dmd = meshpartition(mesh.p,mesh.t,mesh.f,t2t,mesh.tprd,app.elemtype,app.boundaryconditions,mesh.boundaryexpr,mesh.periodicexpr,app.porder,mpiprocs,app.metis);
dmd = meshpartition2(mesh.tprd,ftmp,t2t,app.boundaryconditions,app.nd,app.elemtype,app.porder,mpiprocs,[]);

for i = 1:mpiprocs            
    % create DG nodes
%     if isfield(mesh, 'dgnodes')
    xdg = mesh.dgnodes(:,:,dmd{i}.elempart);
%     else
%         xdg = createdgnodes(mesh.p,mesh.t(:,dmd{i}.elempart),mesh.f(:,dmd{i}.elempart),mesh.curvedboundary,mesh.curvedboundaryexpr,app.porder);    
%         mesh.dgnodes = xdg;
%     end    
    %xdg = createdgnodes(mesh.p,mesh.t(:,dmd{i}.elempart),mesh.f(:,dmd{i}.elempart),mesh.curvedboundary,mesh.curvedboundaryexpr,app.porder);    
    
    [~,cgelcon,rowent2elem,colent2elem,cgent2dgent] = mkcgent2dgent(xdg,1e-6);
%     disp(['Writing initial solution into file ' num2str(i) '...']);
%     if mpiprocs>1
%         fileID1 = fopen(filename + "sol" + string(i) + ".bin",'w');
%     else
%         fileID1 = fopen(filename + "sol" + ".bin",'w');
%     end
%     ndims = zeros(12,1);           
%     ndims(1) = length(dmd{i}.elempart); % number of elements
%     ndims(2) = sum(dmd{i}.facepartpts); % number of faces
%     ndims(3) = size(master.perm,2); % number of faces per element          
%     ndims(4) = master.npe;
%     ndims(5) = master.npf;            
%     ndims(6) = app.nc;
%     ndims(7) = app.ncu;
%     ndims(8) = app.ncq;
%     ndims(9) = app.ncw;
%     ndims(10) = app.nco;
%     ndims(11) = app.nch;
%     ndims(12) = app.ncx;
%     app.nce=10;
%     ndims(13) = app.nce;
%     disp(master.perm);
% 
%     nsize = zeros(20,1);
%     nsize(1) = length(ndims(:));
%     nsize(2) = length(xdg(:));
    %nsize(3) = length(udg(:)); 
%     nsize(4) = length(odg(:));    
% %     nsize(5) = length(wdg(:));    
%     if isfield(mesh, 'udg')        
%         nsize(3) = numel(mesh.udg(:,:,dmd{i}.elempart));
%     end
%     if isfield(mesh, 'vdg')        
%         nsize(4) = numel(mesh.vdg(:,:,dmd{i}.elempart));
%     end
%     if isfield(mesh, 'wdg')        
%         nsize(5) = numel(mesh.wdg(:,:,dmd{i}.elempart));
%     end
%     disp(nsize)
% 
%     fwrite(fileID1,length(nsize(:)),'double',endian);
%     fwrite(fileID1,nsize(:),'double',endian);
%     fwrite(fileID1,ndims(:),'double',endian);
%     fwrite(fileID1,xdg(:),'double',endian);
    %fwrite(fileID1,udg(:),'double',endian);        
%     fwrite(fileID1,odg(:),'double',endian);
%     fwrite(fileID1,wdg(:),'double',endian);
%     if isfield(mesh, 'udg')        
%         fwrite(fileID1,mesh.udg(:,:,dmd{i}.elempart),'double',endian);                
%     end
%     if isfield(mesh, 'vdg')        
%         fwrite(fileID1,mesh.vdg(:,:,dmd{i}.elempart),'double',endian);                
%     end
%     if isfield(mesh, 'wdg')        
%         fwrite(fileID1,mesh.wdg(:,:,dmd{i}.elempart),'double',endian);                
%     end
%     fclose(fileID1);         

    % divide elements and faces into blocks
%     if mpiprocs==1
    ne = length(dmd{i}.elempart);
    [eblks,nbe] = exa_mkelemblocks(ne,2000);
    eblks(3,:) = 0;        
    mf = cumsum([0 dmd{i}.facepartpts]);    
    [fblks,nbf] = exa_mkfaceblocks(mf,dmd{i}.facepartbnd);       
    neb = max(eblks(2,:)-eblks(1,:))+1;
    nfb = max(fblks(2,:)-fblks(1,:))+1;
    app.eblks = eblks;
    app.nbe = nbe;
    app.fblks = fblks;
    app.nbf = nbf;
%     else
%         me = cumsum([0 dmd{i}.elempartpts(1) dmd{i}.elempartpts(2) dmd{i}.elempartpts(3)]);
%         [eblks,nbe] = mkfaceblocks(me,[0 1 2],app.neb);          
%         mf = cumsum([0 dmd{i}.facepartpts]);                 
%         [fblks,nbf] = mkfaceblocks(mf,dmd{i}.facepartbnd,app.nfb);        
%         neb = max(eblks(2,:)-eblks(1,:))+1;
%         nfb = max(fblks(2,:)-fblks(1,:))+1;        
%     end        

    npe = master.npe;
    nfe = size(master.perm,2); 
    facecon1 = reshape(dmd{i}.facecon(:,1,:),[size(dmd{i}.facecon,1) size(dmd{i}.facecon,3)]);
    facecon2 = reshape(dmd{i}.facecon(:,2,:),[size(dmd{i}.facecon,1) size(dmd{i}.facecon,3)]);      
    ind = [];        
    for ii = 1:size(fblks,2)
        if fblks(3,ii)>0
            ind = [ind fblks(1,ii):fblks(2,ii)];
        end
    end         
    
    facecon2(:,ind)=[];        
    [mesh.rowe2f1,mesh.cole2f1,mesh.ent2ind1] = mkdge2dgf(facecon1,npe*length(dmd{i}.elempart));                
    [mesh.rowe2f2,mesh.cole2f2,mesh.ent2ind2] = mkdge2dgf(facecon2,npe*length(dmd{i}.elempart));
    mesh.cole2f1 = mesh.cole2f1 - 1;
    mesh.cole2f2 = mesh.cole2f2 - 1;
    mesh.ent2ind1 = mesh.ent2ind1 - 1;
    mesh.ent2ind2 = mesh.ent2ind2 - 1;
end

% app = writeapp(app,fileapp,endian);                
% 
% mesh.telem = master.telem;
% mesh.tface = master.telem;
% mesh.xpe = master.xpe;
% mesh.xpf = master.xpf;
% for i = 1:mpiprocs   
%     dmd{i}.elem2cpu = [];
%     dmd{i}.elemrecv = [];
%     dmd{i}.nbsd = [];
%     dmd{i}.elemsend = [];
%     dmd{i}.elemsendpts = [];
%     dmd{i}.elemrecvpts = [];
%     dmd{i}.facepartpts = [];
%     dmd{i}.facepartbnd = [];
%     dmd{i}.facecon = [];
% end

