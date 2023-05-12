function A = refinementmatrix(mesh,nref)

nd = mesh.nd;
porder=mesh.porder;
plocal=mesh.plocal;
tlocal=mesh.tlocal;

if isempty(nref), nref=ceil(log2(max(porder,1))); end
if mesh.elemtype==0  
    A0=koornwinder(plocal(:,1:nd),porder);
    plocal = uniref(plocal,tlocal,nref);    
    A=koornwinder(plocal(:,1:nd),porder)/A0;
else
    A0=tensorproduct(plocal(:,1:nd),porder);
    m = porder*(nref+1)+1;     
    if nd==2
        plocal =squaremesh(m,m,1);
    else
        plocal = cubemesh(m,m,m,1);
    end
    A=tensorproduct(plocal(:,1:nd),porder)/A0;  
end

