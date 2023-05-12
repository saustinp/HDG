load result2.mat

a = minf2/min(minf2);
i = find(a>5);
N = i(end);
n = (1:N)';
b1 = a(n);
b2 = lambda2(n);
b3 = kappa2(n);
b4 = b2.*b3;

clear p1 r1;
for i = 1:length(UDG2)
    [x1,y1,p1(:,i)] = getfieldatbou(mesh2,eulereval(UDG2{i},'p',1.4,7));
    [x1,y1,r1(:,i)] = getfieldatbou(mesh2,eulereval(UDG2{i},'r',1.4,7));    
end
maxp = max(p1,[],1);
maxr = max(r1,[],1);
b5 = maxp(n)';
b6 = maxr(n)';
b7 = 0*b6;
b8 = 0*b6;
b9 = 0*b6;
b10 = 0*b6;
for i = 1:length(UDG2)
  dens = eulereval(UDG2{i},'r',1.4,7);
  mach = eulereval(UDG2{i},'M',1.4,7);
  pres = eulereval(UDG2{i},'p',1.4,7);
  b7(i) = max(dens(:));
  b8(i) = max(pres(:));  
  b9(i) = max(mach(:));
  b10(i) = min(pres(:));  
  b6(i) = calerror(ACG2{i},mesh2,master);  
  gam = 1.4;
  Minf = 3.0;
  pinf = 1/(gam*Minf^2);
  ui = 0.5+pinf/(gam-1);
  Hex = gam*ui - (gam-1)*0.5;
  H = eulereval(UDG2{i},'H',1.4,7);
  b5(i) = calerror(H,mesh2,master,@enthalpy);
  b4(i) = max(H(:))-Hex;
end

%tab2 = [n b2 b3 b1 b6 b5 b9 b7 b8 b10];
tab2 = [n b2 b3 b1 b6 b5 b9 b7 b8];



