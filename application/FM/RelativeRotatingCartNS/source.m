function [sr,sr_udg] = source(p,udg,param,time)

nc = size(udg,2);

if nc == 12
    [sr,sr_udg] = source_ns2d(p,udg,param,time);
elseif nc == 20
    [sr,sr_udg] = source_ns3d(p,udg,param,time);
end


