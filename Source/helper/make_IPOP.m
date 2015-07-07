%% make_IPOP.m
% [IPOP] = make_IPOP(basic_map_size)
% generates an initial population size based on making the same initial
% density as for default settings (300 pop/ 2025 area = 0.1481 pop/area)
%
function [IPOP] = make_IPOP(basic_map_size,ipop,),  
[area,~,~] = landscape_measures(basic_map_size);
IDEN = ipop/2025;
IPOP = ceil(IDEN*area);
end