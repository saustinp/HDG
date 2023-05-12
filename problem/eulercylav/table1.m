%load result1.mat

a = minf1/min(minf1);
i = find(a>5);
N = i(end);
n = (1:N)';
b1 = a(n);
b2 = lambda(n);
b3 = kappa(n);
b4 = b2.*b3;

% clear p1 r1;
% for i = 1:length(UDG1)
%     [x1,y1,p1(:,i)] = getfieldatbou(mesh,eulereval(UDG1{i},'p',1.4,7));
%     [x1,y1,r1(:,i)] = getfieldatbou(mesh,eulereval(UDG1{i},'r',1.4,7));    
% end
maxp = max(p1,[],1);
maxr = max(r1,[],1);
b5 = maxp(n)';
b6 = maxr(n)';
b7 = 0*b6;
b8 = 0*b6;
b9 = 0*b6;
b10 = 0*b6;
for i = 1:length(UDG1)
  dens = eulereval(UDG1{i},'r',1.4,7);
  mach = eulereval(UDG1{i},'M',1.4,7);
  pres = eulereval(UDG1{i},'p',1.4,7);
  b7(i) = max(dens(:));
  b8(i) = max(pres(:));  
  b9(i) = max(mach(:));
  b10(i) = min(pres(:));  
  b6(i) = calerror(ACG1{i},mesh,master);
  b5(i) = calerror(eulereval(UDG1{i},'H',1.4,7),mesh,master,@enthalpy);
end

%tab1 = [n b2 b3 b1 b6 b9 b7 b8 b10 b5];
tab1 = [n b2 b3 b1 b6 b5 b9 b7 b8];


