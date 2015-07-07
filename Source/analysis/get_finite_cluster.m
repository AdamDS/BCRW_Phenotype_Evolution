%% get_finite_cluster.m
% Checkes which clusters are finite. Returns a list of values
% answering if a cluster is finite (1 = finite, 0 = infinite along one
% axis, -1 = infinite along both axes).
% -ADS 10*4*12
function [finite_cluster] = get_finite_cluster(tc_of_gen,coords_of_gen,INFRAT)
global SIMOPTS;
N = length(unique(tc_of_gen)); %number of clusters to check
Lx = (((SIMOPTS.basic_map_size(1)*2)-1)*2)-1; 
Ly = (((SIMOPTS.basic_map_size(1)*2)-1)*2)-1;
finite_cluster = ones(N,1); %assume all clusters finite initially
for i = 1:N
  orgs = find(i==tc_of_gen);  %the organisms of current cluster
  x = coords_of_gen(orgs,1);  %x coords
  Mx = max(x);  mx = min(x);  %get min and max x coords
  if ((Mx-mx)/Lx)>=INFRAT
    finite_cluster(i) = finite_cluster(i)-1;  %cluster is infinite along x
  end
  y = coords_of_gen(orgs,2);  %y coords
  My = max(y);  my = min(y);
  if ((My-mx)/Ly)>=INFRAT
    finite_cluster(i) = finite_cluster(i)-1;  %cluster is infinite along y
  end
end
end