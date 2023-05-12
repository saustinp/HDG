function [u,q,uhat,v,iter] = hdg_ma6(master, mesh, app, u, uhat)
% Extend hdg_ma3 to handle linear and nonlinear neumman boundary conditions

[M, Minv, MiB, MiC, MiD, MiE, F, H, G, L, S, D, E, P] = hdg_compute(master, mesh, app);
 
% compute the index mapping
elcon = mesh.elcon;

q = getq(MiB, MiC, u, uhat, elcon);
v = getv(MiD, MiE, u, uhat, elcon);

iter = 0;
itermax = 30;
while (iter<itermax)
  iter = iter + 1;
  [AE, FE, Ua, Wa] = assemble(master, mesh, app, M, MiB, MiC, MiD, MiE, H, G, L, S, v, q, u, uhat);  
  [Ru] = uResidual(master, mesh, app, H, G, L, S, v, q, u, uhat);
  normRu = norm(Ru);
  
  [c, uc] = uaverage(Wa, Ua, P, mesh.elcon);  
  
  [duhat, normR] = linearsystem(AE, FE, elcon, app.neumman, uc, c);
  du = getu(Wa, Ua, duhat, elcon);

  alpha = 1;
  uhat0 = uhat; u0 = u;  
  uhat = uhat0 + alpha*duhat;
  u = u0 + alpha*du;
  q = getq(MiB, MiC, u, uhat, elcon);
  v = getv(MiD, MiE, u, uhat, elcon);  
  [Ru] = uResidual(master, mesh, app, H, G, L, S, v, q, u, uhat);
  normRuNew = norm(Ru);
      
  while normRuNew > normRu 
    alpha = alpha/2;
    uhat = uhat0 + alpha*duhat;
    u = u0 + alpha*du;
    q = getq(MiB, MiC, u, uhat, elcon);
    v = getv(MiD, MiE, u, uhat, elcon);  
    [Ru] = uResidual(master, mesh, app, H, G, L, S, v, q, u, uhat);
    normRuNew = norm(Ru);   
    if alpha < 0.5
      break;
    end
  end  

  [iter alpha normR normRu normRuNew]
  if normR<1e-6
    break;
  end
end

% figure(2); clf; meshplot(mesh,1);    
% mesh2 = mesh; 
% mesh2.dgnodes = q;
% figure(3); clf; meshplot(mesh2,1);

end

function [AE, FE, Ua, Wa, Ru] = assemble(master, mesh, app, M, MiB, MiC, MiD, MiE, H, G, L, S, v, q, u, uhat)

npv = size(M,1);
[npfe, nt] = size(S);
nd = size(MiB,3);
ngv = master.ngv;
npf = master.npf;
ngf = master.ngf;
nfe = size(master.perm,2);

vg = zeros(ngv, nd, nd);
Qa = zeros(npv*nd,npfe);
Ra = zeros(npv*nd,1);
DE = zeros(npv, npv, nd, nd);
DQ = zeros(npv, npv, nd);

Ua = zeros(npv,npfe,nt);
Wa = zeros(npv,nt);
AE = zeros(npfe,npfe,nt);
FE = zeros(npfe,nt);
Ru = zeros(npv,nt);

fbou   = str2func(app.fbou);
source   = str2func(app.source);
shap = squeeze(master.shapvl(:,:,1));
shapfc = squeeze(master.shapfc(:,:,1));
for i = 1:nt % loop over each element    
    Mass = M(:,:,i);    
    M_K = H(:,:,i);
    G_K = G(:,:,i);
    I_K = L(:,:,i);
    S_K = S(:,i);
        
    % form F_K
    [xg, ~, jac] = volgeom(master.shapmv,reshape(mesh.dgnodes(:,:,i),[npv 1 nd]));      
    ug = master.shapmv(:,:,1)*u(:,1,i);
    qg = master.shapmv(:,:,1)*q(:,:,i);
    [sg, sg_udg] = source(xg, [ug qg], app.param);            
    %sg = source(xg, [xg xg(:,1)], app.param);    
    for l = 1:nd
      vg(:,:,l) = master.shapmv(:,:,1)*v(:,:,l,i);
    end        
    tg = sqrt(2*sg + vg(:,1,1).^2 + vg(:,2,2).^2 + vg(:,1,2).^2 + vg(:,2,1).^2);
    fg = vg(:,1,1) + vg(:,2,2) - tg;
    F_K = shap*(master.gwvl.*jac.*fg);        
    Ru(:,i) = F_K;
    
    for l = 1:nd
      DQ(:,:,l)  = shap*diag(master.gwvl.*jac.*sg_udg(:,1,1+l)./tg)*shap';    
      for m = 1:nd
        DE(:,:,l,m)  = shap*diag(master.gwvl.*jac.*vg(:,l,m)./tg)*shap';    
      end
    end    
    
    Di = -Mass*MiD(:,:,1,1,i);
    Ei = Mass*MiE(:,:,1,1,i);
    for l = 2:nd
      Di = Di - Mass*MiD(:,:,l,l,i);
      Ei = Ei + Mass*MiE(:,:,l,l,i);
    end        

    Di = Di - (DQ(:,:,1)*MiB(:,:,1,i) + DQ(:,:,2)*MiB(:,:,2,i));
    Ei = Ei - (DQ(:,:,1)*MiC(:,:,1,i) + DQ(:,:,2)*MiC(:,:,2,i));
    
    Di = Di + DE(:,:,1,1)*MiD(:,:,1,1,i) + DE(:,:,2,2)*MiD(:,:,2,2,i)  ...
            + DE(:,:,1,2)*MiD(:,:,1,2,i) + DE(:,:,2,1)*MiD(:,:,2,1,i);
    Ei = Ei - DE(:,:,1,1)*MiE(:,:,1,1,i) - DE(:,:,2,2)*MiE(:,:,2,2,i)  ...
            - DE(:,:,1,2)*MiE(:,:,1,2,i) - DE(:,:,2,1)*MiE(:,:,2,1,i);    
        
    % Di*U = Ei*Uhat + F_K
    Ua(:,:,i) = Di\Ei;    
    Wa(:,i) = Di\F_K;    
    
    imap = mesh.elcon(:,i);        
    Rh = S_K - M_K*uhat(imap) - G_K*reshape(q(:,:,i),[npv*nd 1]) + I_K*u(:,i);
    
    nbc = 0; % check boundary conditions
    for j = 1:nfe 
      ibf = -mesh.f(abs(mesh.t2f(i,j)),end);
      if (ibf>0) && (app.bcm(ibf) > 1)          
        nbc = 1;        
      end
    end

    if nbc == 1
      [xf, nlg, jacf] = facegeom(master.shapmf,reshape(mesh.dgnodes(:,:,i),[npv 1 nd]),master.perm);        
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
            qg = reshape(master.shapmf(:,:,1)*qn,[ngf nd]);    
            [gb, dgb] = fbou(xgf,nl,pars,app.bcm(ibf),qg);                  
            G_K(J,:) = 0;
            for l=1:nd              
                G_K(J,(l-1)*npv+I) = G_K(J,(l-1)*npv+I) + shapfc*diag(dws.*dgb(:,l))*shapfc';             
            end                  
            S_K(J,:) = shapfc*(dws.*(-gb));  
            Rh(J,:) = S_K(J,:) - M_K(J,:)*uhat(imap) + I_K(J,:)*u(:,i);
          end
      end
    end            
        
    for l = 1:nd
      Qa((l-1)*npv+1:l*npv,:) = MiC(:,:,l,i) - MiB(:,:,l,i)*Ua(:,:,i);
      Ra((l-1)*npv+1:l*npv,:) = MiB(:,:,l,i)*Wa(:,i);
    end    
    H_K = M_K + G_K*Qa - I_K*Ua(:,:,i);
    R_K = G_K*Ra + I_K*Wa(:,i);      
    
    AE(:,:,i) = H_K;
    FE(:,i) = R_K + Rh;        
end

end

function [Ru] = uResidual(master, mesh, app, H, G, L, S, v, q, u, uhat)

npv = master.npv;
nt = size(mesh.dgnodes,3);
nd = master.nd;
ngv = master.ngv;
npf = master.npf;
ngf = master.ngf;
npfe = size(S,1);
nfe = size(master.perm,2);

vg = zeros(ngv, nd, nd);
Ru = zeros(npv,nt);
Rh = zeros(npfe,nt);

fbou   = str2func(app.fbou);
source   = str2func(app.source);
shap = squeeze(master.shapvl(:,:,1));
shapfc = squeeze(master.shapfc(:,:,1));
for i = 1:nt % loop over each element            
    % form F_K
    [xg, ~, jac] = volgeom(master.shapmv,reshape(mesh.dgnodes(:,:,i),[npv 1 nd]));      
    %sg = source(xg, [xg xg(:,1)], app.param);    
    ug = master.shapmv(:,:,1)*u(:,1,i);
    qg = master.shapmv(:,:,1)*q(:,:,i);
    sg = source(xg, [ug qg], app.param);            
    for l = 1:nd
      vg(:,:,l) = master.shapmv(:,:,1)*v(:,:,l,i);
    end        
    tg = sqrt(2*sg + vg(:,1,1).^2 + vg(:,2,2).^2 + vg(:,1,2).^2 + vg(:,2,1).^2);
    fg = vg(:,1,1) + vg(:,2,2) - tg;
    F_K = shap*(master.gwvl.*jac.*fg);        
    Ru(:,i) = F_K;    
    
    S_K = S(:,i);
    M_K = H(:,:,i);
    G_K = G(:,:,i);
    I_K = L(:,:,i);    
    imap = mesh.elcon(:,i);        
    Rh(:,i) = S_K - M_K*uhat(imap) - G_K*reshape(q(:,:,i),[npv*nd 1]) + I_K*u(:,i);    
    
    nbc = 0; % check boundary conditions
    for j = 1:nfe 
      ibf = -mesh.f(abs(mesh.t2f(i,j)),end);
      if (ibf>0) && (app.bcm(ibf) > 1)          
        nbc = 1;        
      end
    end

    if nbc == 1
      [xf, nlg, jacf] = facegeom(master.shapmf,reshape(mesh.dgnodes(:,:,i),[npv 1 nd]),master.perm);        
      xf = reshape(xf, [ngf nfe nd]);
      nlg = reshape(nlg, [ngf nfe nd]);
      jacf = reshape(jacf, [ngf nfe]);    
      for j = 1:nfe % loop over each face of the element i
          ibf = -mesh.f(abs(mesh.t2f(i,j)),end);
          if (ibf>0) && (app.bcm(ibf) > 1)          
            J = ((j-1)*npf+1):j*npf;        
            pars = app.bcd(ibf,:);
            
            nl = reshape(nlg(:,j,:),[ngf nd]);
            xgf = reshape(xf(:,j,:),[ngf nd]);
            
            dws = master.gwfc.*jacf(:,j);                             
            qn = reshape(q(master.perm(:,j),:,i),[npf nd]);    
            qg = reshape(master.shapmf(:,:,1)*qn,[ngf nd]);    
            gb = fbou(xgf,nl,pars,app.bcm(ibf),qg);                  
            S_K(J,:) = shapfc*(dws.*(-gb));  
            Rh(J,i) = S_K(J,:) - M_K(J,:)*uhat(imap) + I_K(J,:)*u(:,i);
          end
      end
    end                
end

elcon = mesh.elcon;
il = zeros(npfe,npfe,nt);
jl = zeros(npfe,npfe,nt);
for i=1:nt
    con = repmat((elcon(:,i)'-1),1,1)+ones(1,npfe);
    con = reshape(con,npfe,1);
    il(:,:,i) = repmat(con ,1,npfe);
    jl(:,:,i) = repmat(con',npfe,1);        
end

R = sparse(reshape(il(:,1,:),(npfe)*nt,1),ones((npfe)*nt,1),reshape(Rh,(npfe)*nt,1));                    
Ru = [Ru(:); R(:)];

end

function [uhat, normR] = linearsystem(AE, FE, elcon, neumann, uc, c)

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

normR = norm(R);

if neumann==1
  R(end+1) = -c;
  H(end+1,:) = uc;  
  H(:,end+1) = [uc(:); 0*max(abs(uc(:)))];
  ut = H\R;
  uhat = ut(1:end-1);    
else  
  uhat = H\R;
end


end


