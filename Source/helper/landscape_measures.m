%% landscape_measures.m
% [A,xl,yl] = landscape_measures(basic_map_size)
% A is the area
% xl is the landscape horizontal length
% yl is the landscape vertical length
%
function [A,xl,yl] = landscape_measures(basic_map_size),  
xl = (((2*basic_map_size(1))-1)*2)-1;
yl = (((2*basic_map_size(2))-1)*2)-1;
A = xl*yl;
end