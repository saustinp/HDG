function q0 = initq_func_z(dgnodes, param)
    r = dgnodes(:,1,:);
    z = dgnodes(:,2,:);
    
    q0 = -(-(7547.*exp(17./10 - (5000.*r.^2)./19 - (5000.*(z - 17./25).^2)./19).*((10000.*z)./19 - 6800./19))./10);
end