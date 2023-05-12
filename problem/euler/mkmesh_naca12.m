function mesh = mkmesh_naca12(porder)

[p,t] = naca0012();

elemtype = 0;
nodetype = 1;
bndexpr = {'all(sqrt(sum(p.^2,2))<2)','all(sqrt(sum(p.^2,2))>2)'};  

mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);

fb = 1;
mesh.fcurved = (mesh.f(:,end)==-fb);
ic = mesh.fcurved;
mesh.tcurved = false(size(mesh.t,1),1);
mesh.tcurved(mesh.f(ic,end-1)) = true;
mesh.dgnodes = makedgnodes(mesh,@naca12dist,fb);

