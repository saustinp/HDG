function [u,q,uhat,v,iter] = hdg_mah(master, mesh, app, rho, u, uhat)
% Extend hdg_ma5 to deal with field density

[M, ~, MiB, MiC, MiD, MiE, F, H, G, L, S, ~, ~, P] = hdg_compute(master, mesh, app);
 
% compute the index mapping
elcon = mesh.elcon;

q = getq(MiB, MiC, u, uhat, elcon);
v = getv(MiD, MiE, u, uhat, elcon);
  
[AE, ~, Ua, ~, Da] = assemble(M, MiB, MiC, MiD, MiE, F, H, G, L, S);

mesh2 = mkcgmesh(mesh);
mesh2.ib = [];
mesh2.in = 1:size(mesh2.p2,1);

sz = size(mesh.dgnodes(:,1:master.nd,:));
qg = zeros(master.nge*mesh.ne,master.nd);

vold = v;
qold = q;
itermax = 100;
iter = 0;
while (iter<itermax)
  
  tm = master.shapmv(:,:,1)*reshape(q(:,1,:),sz(1),[]);   
  qg(:,1) = tm(:);
  tm = master.shapmv(:,:,1)*reshape(q(:,2,:),sz(1),[]); 
  qg(:,2) = tm(:);
  if iter==0    
    [elist, xi] = locatexinmesh(mesh, qg, [], 1e-4);
  else
    [elist, xi] = locatexinmesh(mesh, qg, elist, 1e-4);
  end  
  rhoxi = evalfield(mesh, rho, elist, xi);
  rhoxi = reshape(rhoxi, [master.nge mesh.ne]); 
  
  [AE, FE, Wa] = assembleRHS(master, mesh, app, AE, rhoxi, u, q, v, MiB, MiC, Da, Ua, H, G, L, S);
  [c, uc] = uaverage(Wa, Ua, P, mesh.elcon);  
  uhat = linearsystem(AE, FE, elcon, app.neumman, uc, c);  
  u = getu(Wa, Ua, uhat, elcon);  
  q = getq(MiB, MiC, u, uhat, elcon);
  v = getv(MiD, MiE, u, uhat, elcon);
  
%   for d = 1:2
% %    figure(1); clf; scaplot(mesh2,q(:,d,:)); axis equal; axis tight; colormap jet;
%     s = cgpoisson(mesh2, master, q(:,d,:), [1e-4 1.0]);    
%     s = s(mesh2.t2');
%     q(:,d,:) = reshape(s, size(q(:,d,:)));
% %     figure(2); clf; scaplot(mesh2,q(:,d,:)); axis equal; axis tight; colormap jet;  
% %     pause
%     for j = 1:2
%       s = cgpoisson(mesh2, master, v(:,d,j,:), [1e-4 1.0]);    
%       s = s(mesh2.t2');
%       v(:,d,j,:) = reshape(s, size(v(:,d,j,:)));
%     end
%   end
  
  if norm(v(:)-vold(:)) <= 1e-4
    break;
  end
        
  [iter norm(v(:)-vold(:)) norm(q(:)-qold(:))]
  vold = v;
  qold = q;
  iter = iter + 1;
  
  if rem(iter,1)==0 %TODO: CHANGE BACK
    mesh2 = mesh;  mesh2.dgnodes = q;
    figure(2); clf; meshplot(mesh2,1);    
  end
  
%   figure(1); clf; scaplot(mesh,v(:,1,2,:)); axis equal; axis tight; colormap jet;
%   figure(2); clf; scaplot(mesh,v(:,2,1,:)); axis equal; axis tight; colormap jet;  
end

% figure(2); clf; meshplot(mesh,1);    

end


function [AE, FE, Ua, Wa, Da] = assemble(M, MiB, MiC, MiD, MiE, F, H, G, L, S)

npv = size(M,1);
[npfe, nt] = size(S);
nd = size(MiB,3);

Qa = zeros(npv*nd,npfe);
Ra = zeros(npv*nd,1);

Da = zeros(npv,npv,nt);
Ua = zeros(npv,npfe,nt);
Wa = zeros(npv,nt);
AE = zeros(npfe,npfe,nt);
FE = zeros(npfe,nt);

for i = 1:nt % loop over each element    
    Mass = M(:,:,i);    
    F_K = F(:,i);
    M_K = H(:,:,i);
    G_K = G(:,:,i);
    I_K = L(:,:,i);
    S_K = S(:,i);
        
    Di = -Mass*MiD(:,:,1,1,i);
    Ei = Mass*MiE(:,:,1,1,i);
    for l = 2:nd
      Di = Di - Mass*MiD(:,:,l,l,i);
      Ei = Ei + Mass*MiE(:,:,l,l,i);
    end        
    % Di*U = Ei*Uhat + F_K
    Da(:,:,i) = Di;
    Ua(:,:,i) = Di\Ei;    
    Wa(:,i) = Di\F_K;    
    
    for l = 1:nd
      Qa((l-1)*npv+1:l*npv,:) = MiC(:,:,l,i) - MiB(:,:,l,i)*Ua(:,:,i);
      Ra((l-1)*npv+1:l*npv,:) = MiB(:,:,l,i)*Wa(:,i);
    end    
    H_K = M_K + G_K*Qa - I_K*Ua(:,:,i);
    R_K = G_K*Ra + I_K*Wa(:,i);      
    
    AE(:,:,i) = H_K;
    FE(:,i) = R_K + S_K;    
end

end

function [AE, FE, Wa] = assembleRHS(master, mesh, app, AE, rho, u, q, v, MiB, MiC, Da, Ua, H, G, L, S)

npv = size(Da,1);
[npfe, nt] = size(S);
nd = size(MiB,3);
ngv = master.ngv;
npf = master.npf;
ngf = master.ngf;
nfe = size(master.perm,2);

vg = zeros(ngv, nd, nd);
Qa = zeros(npv*nd,npfe);
Ra = zeros(npv*nd,1);
Wa = zeros(npv,nt);
FE = zeros(npfe,nt);

fbou   = str2func(app.fbou);
%source = str2func(app.source);
shapfc = squeeze(master.shapfc(:,:,1));
shap = squeeze(master.shapvl(:,:,1));
for i = 1:nt % loop over each element    
    G_K = G(:,:,i);
    I_K = L(:,:,i);
    S_K = S(:,i);
            
    % form F_K
    [~, ~, jac] = volgeom(master.shapmv,reshape(mesh.dgnodes(:,1:nd,i),[npv 1 nd]));      
    %ug = master.shapmv(:,:,1)*u(:,1,i);
    %qg = master.shapmv(:,:,1)*q(:,:,i);
    %sg = source(xg, [ug qg], app.param);        
    sg = app.param(1)./rho(:,i);
    for l = 1:nd
      vg(:,:,l) = master.shapmv(:,:,1)*v(:,:,l,i);
    end            
    fg = -sqrt(2*sg + vg(:,1,1).^2 + vg(:,2,2).^2 + vg(:,1,2).^2 + vg(:,2,1).^2);    
    F_K = shap*(master.gwvl.*jac.*fg);    
        
    Wa(:,i) = Da(:,:,i)\F_K;        
    for l = 1:nd
      Ra((l-1)*npv+1:l*npv,:) = MiB(:,:,l,i)*Wa(:,i);
    end    
    
    nbc = 0; % check boundary conditions
    for j = 1:nfe 
      ibf = -mesh.f(abs(mesh.t2f(i,j)),end);
      if (ibf>0) && (app.bcm(ibf) > 1)          
        nbc = 1;        
      end
    end
    
    if nbc == 1
      M_K = H(:,:,i);      
      [xf, nlg, jacf] = facegeom(master.shapmf,reshape(mesh.dgnodes(:,1:nd,i),[npv 1 nd]),master.perm);        
      xf = reshape(xf, [ngf nfe nd]);
      nlg = reshape(nlg, [ngf nfe nd]);
      jacf = reshape(jacf, [ngf nfe]);    
      for j = 1:nfe % loop over each face of the element i
          ibf = -mesh.f(abs(mesh.t2f(i,j)),end);
          if (ibf>0) && (app.bcm(ibf) > 1)          
            I = master.perm(:,j,1);
            J = ((j-1)*npf+1):j*npf;        
            pars = app.bcd(ibf,:);
            nl = reshape(nlg(:,j,:),[ngf nd]);
            xgf = reshape(xf(:,j,:),[ngf nd]);
            
            dws = master.gwfc.*jacf(:,j);                             
            qn = reshape(q(master.perm(:,j),:,i),[npf nd]);    
            %qn(:,1).^2 + qn(:,2).^2
            qg = reshape(master.shapmf(:,:,1)*qn,[ngf nd]);                   
            G_K(J,:) = 0;
%             for l=1:nd              
%                 tmx = shapfc*diag(dws.*qg(:,l)/pars(l))*shapfc';            
%                 G_K(J,(l-1)*npv+I) = G_K(J,(l-1)*npv+I) - 2*tmx;   
%             end                  
            %gb = fbou(xgf,nl,app.param,app.bcm(ibf));   
            %S_K(J,:) = shapfc*(dws.*(gb - qg(:,1).^2/pars(1) - qg(:,2).^2/pars(2)));             
            [gb, dgb] = fbou(xgf,nl,pars,app.bcm(ibf),qg);               
            for l=1:nd                               
                G_K(J,(l-1)*npv+I) = G_K(J,(l-1)*npv+I) + shapfc*diag(dws.*dgb(:,l))*shapfc';             
            end                                          
            S_K(J,:) = shapfc*(dws.*(-gb + dgb(:,1).*qg(:,1) + dgb(:,2).*qg(:,2)));             
          end
      end
      for l=1:nd      
        Qa((l-1)*npv+1:l*npv,:) = MiC(:,:,l,i) - MiB(:,:,l,i)*Ua(:,:,i);
      end
      AE(:,:,i) = M_K + G_K*Qa - I_K*Ua(:,:,i);           
    end        
    
    FE(:,i) = G_K*Ra + I_K*Wa(:,i) + S_K;    
end

end


function uhat = linearsystem(AE, FE, elcon, neumann, uc, c)

[npfe, nt] = size(FE);

il = zeros(npfe,npfe,nt);
jl = zeros(npfe,npfe,nt);
for i=1:nt
    con = repmat((elcon(:,i)'-1),1,1)+ones(1,npfe);
    con = reshape(con,npfe,1);
    il(:,:,i) = repmat(con ,1,npfe);
    jl(:,:,i) = repmat(con',npfe,1);        
end
H = sparse(reshape(il,(npfe)^2*nt,1),reshape(jl,(npfe)^2*nt,1),reshape(AE,(npfe)^2*nt,1));
R = sparse(reshape(il(:,1,:),(npfe)*nt,1),ones((npfe)*nt,1),reshape(FE,(npfe)*nt,1));                    

if neumann==1
  R(end+1) = -c;
  H(end+1,:) = uc;  
  H(:,end+1) = [uc(:); 4*max(abs(uc(:)))];
  ut = H\R;
  uhat = ut(1:end-1);      
else  
  uhat = H\R;
end

end
