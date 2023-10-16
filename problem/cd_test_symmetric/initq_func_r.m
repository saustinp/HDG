function q0 = initq_func_r(dgnodes, param)
    r = dgnodes(:,1,:);
    z = dgnodes(:,2,:);

    q0 = -(-(7547000.*r.*exp(17./10 - (5000*r.^2)./19 - (5000.*(z - 17./25).^2)./19))./19);
end