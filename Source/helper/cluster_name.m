%% cluster_name.m
% [clus_name] = cluster_name(base_name)
% Corrects base_name when limit isn't default (==3)
%
% -ADS 5*3*13
function [clus_name] = cluster_name(base_name), 
global SIMOPTS;
if SIMOPTS.limit~=3,  clus_name = [base_name int2str(SIMOPTS.limit) '_limit_'];
else, clus_name = base_name;  end
end