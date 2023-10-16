function u0 = initu_func(dgnodes, param)
    r = dgnodes(:,1,:);
    z = dgnodes(:,2,:);

    u0 = 754.7 *exp(((z-.68).^2 + r.^2)./-0.0038 + 1.7);
end