function [UDG,UH]=hdg_solve(master,mesh,app,UDG,UH,SH)
%HDG_SOLVE Solve using the HDG method and Newton iteration
%   [UH,QH,PH,UHAT] = HDG_SOLVE(MASTER,MESH,UDG,UH,SH,APP)
%
%      MASTER:                  Master structure
%      MESH:                    Mesh structure
%      UH(NPL,NC,NE):           Vector of unknowns (initial guess)
%      QH(NPL,NC,2,NE):         Vector of gradients of U (initial guess)
%      PH(NPL,1,NE):            Vector of pressure (initial guess)
%      UHAT(NC,3*NPS,NF):       Vector of U_hat's (initial guess)
%      APP:                     Application structure
%      SH:                      Source term independent of U (for time dependent mode only)
%
%      UH(NPL,NC,NE):           Vector of unknowns
%      QH(NPL,NC,2,NE):         Vector of gradients of U
%      PH(NPL,1,NE):            Vector of pressure
%      UHAT(NCH,3*NPS,NF):      Vector of U_hat's
%
%      NPL:                     Number of DG nodes within an element
%      NC:                      Number of conservation equations solved (components)
%      NPS:                     Number of HDG nodes per edge (porder+1)
%      NE:                      Number of elements
%      NF:                      Number of faces
%


npf = master.npf;
npv = master.npv;
ne  = mesh.ne;
%nf  = mesh.nf;
nfe = size(master.perm,2);
nc  = app.nc;
ncu = app.ncu;
%nd  = app.nd;
nsiz = mesh.nsiz;

if isfield(app,'fbou') == 0
    app.fbou = 'fbou';
end
if isfield(app,'source') == 0
    app.source = 'source';
end
if isfield(app,'denseblock') == 0
    app.denseblock = 0;
end
if isfield(app,'linear') == 0
    app.linear = 0;
end
if isfield(app,'adjoint') == 0
    app.adjoint = 0;
end
if isfield(app,'getdqdg') == 0
    app.getdqdg = 1;
end

if app.getdqdg == 0
    nco = ncu;
else
    nco = nc;
end
if isempty(SH)
    SH = zeros(npv,nc,ne);
end
if app.wave==0
    SH(:,ncu+1:end,:)=0;
else
    if app.flg_p
        SH(:,ncu+1,:)=0;
    end
end

if app.denseblock ~= 0
    f2f = mkf2f(mesh.f, mesh.t2f);
end

if min(master.perm(:))==0
    master.perm = master.perm + 1;
end
if min(master.permgeom(:))==0
    master.permgeom = master.permgeom + 1;
end

tic
[K,F,DUDG,DUDG_DUH] = hdg_assemble(master,mesh,app,UDG,UH,SH);  
toc

if app.adjoint == 1 
    if app.denseblock == 0
        UHE  = full(reshape(K\F,ncu,nsiz));        
    else
        UHE = gmres(K, F, f2f);
        UHE = full(reshape(UHE,ncu,nsiz));    
    end
    if isempty(DUDG) == 0
        UHE = reshape(full(UHE(:,mesh.elcon(:))),ncu,nfe*npf,ne);    
        UDG = DUDG + reshape(mapContractK(DUDG_DUH,UHE,[1 2],[3 4],5,[1 2],[],3),[npv nco ne]);
    else        
        UDG = hdg_local(master,mesh,app,UDG,UH,SH,UHE);        
    end
    UH = UHE;
    return;
end

it   = 0;
duh  = 1e6;
NewtonTol = 1e-8;
while duh > NewtonTol && it < 16
                
    if app.denseblock==0
        DUH = full(reshape(K\F,ncu,nsiz));   
%         DUH(:,1:3*npf)
%         pause
    else
        DUH = gmres(K, F, f2f);
        DUH = full(reshape(DUH,ncu,nsiz));    
    end    
                
    if isempty(DUDG)==0
        DUHE = reshape(full(DUH(:,mesh.elcon(:))),ncu,nfe*npf,ne);    
        DUDGT = DUDG + reshape(mapContractK(DUDG_DUH,DUHE,[1 2],[3 4],5,[1 2],[],3),[npv nco ne]);
    else % compute DUDGT by solving the local problems
        DUDGT = hdg_local(master,mesh,app,UDG,UH,SH,DUH);
    end
    
    if app.linear == 1
        UDG(:,1:nco,:)  = UDG(:,1:nco,:) + DUDGT;                     
        UH   = UH + DUH;   
        if app.getdqdg==0
            QDG = getq(master, mesh, UDG, UH, SH, app.fc_q);
            UDG(:,ncu+1:end,:)=QDG;
        end        
        return; 
    end
    
    it = it + 1;    
    fprintf('Newton iteration :  %d\n', it);
    
    % Find stepsize for damped Newton    
    UDG0 = UDG;
    UH0  = UH;    
    duh0 = norm(F(:));   
    alfa = 1;
    while 1        
        UDG(:,1:nco,:)  = UDG0(:,1:nco,:) + alfa*DUDGT;                     
        UH   = UH0 + alfa*DUH;   
        
        if app.getdqdg==0
            QDG = getq(master, mesh, UDG, UH, SH, app.fc_q);
            UDG(:,ncu+1:end,:)=QDG;
        end        
        
        [K,F,DUDG,DUDG_DUH] = hdg_assemble(master,mesh,app,UDG,UH,SH);        
        duh  = norm(F(:));                 
        
        if duh>duh0
            alfa=alfa/2;
            fprintf(' alpha = %f\n',alfa);
            if alfa<1e-2, break; end
        else
            break;
        end
    end
       
    %save tmp.mat UDG UH;
    
%     alpha = 100; beta = 0.1; href = 0.12; hk = 0.001;
%     div = divergence(UDG, href);
%     s = ucg2udg(udg2ucg(div, mesh.cgent2dgent, mesh.rowent2elem), mesh.cgelcon);
%     s = cgpoisson(mesh, master, s, [hk 1.0]);        
%     a = (s-beta).*(atan(alpha*(s-beta))/pi + 0.5) - atan(alpha)/pi + 0.5;    
%     a = reshape(a(mesh.t2'), master.npe, 1, mesh.ne);
%     mesh.dgnodes(:,mesh.nd+1,:) = mesh.av*a;
        
    fprintf('Old residual: %e,   New residual: %e    %e\n', [duh0 duh alfa]);       
    %figure(2); clf; scaplot(mesh,eulereval(UDG(:,1:app.ncu,:),'M',1.4),[],1,1); axis off; %axis equal; axis tight;
    
%     e0=8.8540e-12;E0=60000;Etstar=1;L0=1;
%     Et = Etstar;     
%     Ex = UDG(:,3,:); Ey = UDG(:,5,:)+Et; E = sqrt(Ex.^2+Ey.^2);
%     [dgnodes,u] = potentialfield(master,mesh,E,-5,1);    
%     [~,rho] = potentialfield(master,mesh,UDG(:,2,:),-5,1);    
%     x = dgnodes(:,1,:);x=x(:);
%     y = dgnodes(:,2,:);y=y(:);
%     f = u(:,1,:);f=f(:);
%     q = rho(:,1,:);q=q(:);
%     t = cart2pol(x,y);
%     i = find(t<0);
%     t(i) = 2*pi+t(i);
%     [tj,jj] = sort(t);
%     figure(1);clf;
%     plot(tj,E0*f(jj)/1e6,'b-','LineWidth',1.5);        
%     figure(2);clf;
%     plot(tj,1e6*e0*E0/L0*q(jj),'-b','LineWidth',1.5);    
%     
%     figure(7); clf; scaplot(mesh,e0*E0/L0*UDG(:,2,:)*1e6,[],2); 
%     xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
%     colorbar('FontSize',15); set(gca,'FontSize',18); box off;
%     axis equal; axis tight; colormap jet; 
%     figure(6); clf; scaplot(mesh,(L0*E0)*UDG(:,1,:)/1e6,[],2); 
%     xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
%     colorbar('FontSize',15); set(gca,'FontSize',18); box off;
%     axis equal; axis tight; colormap jet;     
%     pause(1);
end
