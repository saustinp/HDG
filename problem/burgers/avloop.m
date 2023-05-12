function [UDG, UH, ACG, mine] = avloop(master, mesh, app, UDG0, av, href, hk, alpha)

nd = mesh.nd;
ncu = app.ncu;
porder = mesh.porder;

if nd == 1
    np1 = porder;
elseif nd==2
    if mesh.elemtype==0
        np1 = porder*(porder+1)/2;
    else
        np1 = porder*porder;
    end
else
    if mesh.elemtype==0
        np1 = porder*(porder+1)*(porder+2)/6;
    else
        np1 = porder*porder*porder;
    end
end

UDG = cell(length(av),1);
UH = cell(length(av),1);
ACG = cell(length(av),1);
UDG{1} = UDG0;
UH{1} = inituhat(master,mesh.elcon,UDG{1},ncu);

% HDG solver for constant viscosity field
mesh.dgnodes(:,nd+1,:) = av(1);
[UDG{1}, UH{1}] = hdg_solve(master,mesh,app,UDG{1},UH{1},[]);

A = tensorproduct(master.plocvl, mesh.porder);
[~,cgelcon,rowent2elem,~,cgent2dgent] = mkcgent2dgent(mesh.dgnodes(:,1:nd,:),1e-8);

mine = zeros(length(av),1);
u = A\squeeze(UDG{1}(:,1,:));
u1 = u(1:np1,:);
e1 = sum(abs(u1),1);
e = sum(abs(u),1);
idx = e>1e-5;
mine(1) = min(e1(idx)./e(idx));

for i = 2:length(av)
    div = divergence(UDG{i-1}, href);
    s = ucg2udg(udg2ucg(div, cgent2dgent, rowent2elem), cgelcon);
    s = cgpoisson(mesh, master, s, [hk 1.0]);        
    a = s.*(atan(alpha*s)/pi + 0.5) - atan(alpha)/pi + 0.5;
    mesh.dgnodes(:,nd+1,:) = av(i)*reshape(a(mesh.t2'), master.npe, 1, mesh.ne);                
    ACG{i-1} = mesh.dgnodes(:,nd+1,:);               
    
    [UDG{i}, UH{i}] = hdg_solve(master,mesh,app,UDG{i-1},UH{i-1},[]);
            
    u = A\squeeze(UDG{i}(:,1,:));
    u1 = u(1:np1,:);
    e1 = sum(abs(u1),1);
    e = sum(abs(u),1);
    idx = e>1e-5;
    mine(i) = min(e1(idx)./e(idx));
    
    [i av(i) mine(i)]       
end

div = divergence(UDG{i}, href);
s = ucg2udg(udg2ucg(div, cgent2dgent, rowent2elem), cgelcon);
a = s.*(atan(alpha*s)/pi + 0.5) - atan(alpha)/pi + 0.5;
a = cgpoisson(mesh, master, a, [hk 1.0]);        
ACG{i} = av(i)*reshape(a(mesh.t2'), master.npe, 1, mesh.ne);        

