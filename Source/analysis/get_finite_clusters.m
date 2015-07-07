%% get_finite_clusters.m
% function [finite_clusters] = get_finite_clusters(population,num_clusters,trace_cluster,...
%                                                  trace_x,trace_y,INFRAT)
% Checkes which clusters are finite. Returns a list of values
% answering if a cluster is finite (1 = finite, 0 = infinite along one
% axis, -1 = infinite along both axes).
% -ADS 3*5*12
function [finite_clusters] = get_finite_clusters(p,nc,tc,tx,ty,INFRAT),
global SIMOPTS;
ngen = length(p(p>=SIMOPTS.limit)); %number of generations
Lx = (((SIMOPTS.basic_map_size(1)*2)-1)*2)-1; 
Ly = (((SIMOPTS.basic_map_size(1)*2)-1)*2)-1;
finite_clusters = ones(sum(nc),1); %assume all clusters finite initially
i = 0;
for gen = 1:ngen, 
  v = sum(p(1:gen)); u = v -p(gen) +1;
  tcg = tc(u:v);
  xg = tx(u:v); 
  yg = ty(u:v);
  ncg = nc(gen);
  cv = sum(nc(1:gen));  cu = cv -nc(gen) +1;
  for c = 1:ncg, 
    i = i +1;
    orgs = find(c==tcg);  %the organisms of current cluster
    x = xg(orgs);  %x coords
    Mx = max(x);  mx = min(x);  %get min and max x coords
    if ((Mx-mx)/Lx)>=INFRAT
      finite_clusters(i) = finite_clusters(i)-1;  %cluster is infinite along x
    end
    y = yg(orgs);  %y coords
    My = max(y);  my = min(y);
    if ((My-my)/Ly)>=INFRAT
      finite_clusters(i) = finite_clusters(i)-1;  %cluster is infinite along y
    end
  end %i
end %gen
end