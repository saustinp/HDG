function a = avfield(mesh, master, UDG)

href = 1;
alpha = 1000;    
div = divergence(UDG, href);
divmax = max(div(:));
s = limiting(div,0,divmax/4,alpha,0);
%s = ucg2udg(udg2ucg(s, mesh.cgent2dgent, mesh.rowent2elem), mesh.cgelcon);

s = cgpoisson(mesh, master, s, [0.0001 1.0]);    
s = s(mesh.t2');
s = reshape(s, size(UDG(:,1,:)));
a = s/max(s(:));

%rho = 1 + 20*a;
figure(1); clf; 
scaplot(mesh, a(:,1,:),[],2,1); axis on; axis equal; axis tight;

