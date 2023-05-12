function uh = cuda_getuhat(ib, uh, udg1, udg2, uinf, xdg1, nlg1, param, time, ng, ncu, nd, nc)

uh = reshape(uh,[ng,ncu]);
if ib==0
    for i=1:ncu
        uh(:,i) = (udg1(:,i) + udg2(:,i))/2 + sign(nlg1(:,1)+nlg1(:,2)).*(udg1(:,i) - udg2(:,i))/2;
    end
elseif ib>0
    uh(:,1:ncu) = ldgubou(ib,uinf,nlg1,xdg1,udg1,param,time);
end
