%% get_cluster_spans.m
% function [cluster_spans] = get_cluster_spans(population,num_clusters,trace_cluster,..
%                                                  trace_x,trace_y)
% Determines the landscape spanning percentage of each cluster along both [x,y] coordinates.
% -ADS 9*25*13
function [cluster_spans] = get_cluster_spans(p,nc,tc,tx,ty),
global SIMOPTS;
ngen = length(p(p>=SIMOPTS.limit)); %number of generations
Lx = (((SIMOPTS.basic_map_size(1)*2)-1)*2)-1; 
Ly = (((SIMOPTS.basic_map_size(1)*2)-1)*2)-1;
cluster_spans = zeros(sum(nc),2); %
i = 0;
for gen = 1:ngen, 
  [u,v] = gen_index(p,gen);
  [cu,cv] = gen_index(nc,gen);
  tcg = tc(u:v);
  xg = tx(u:v); 
  yg = ty(u:v);
  ncg = nc(gen);
  for c = 1:ncg, 
    C = c+cu-1;
    i = i +1;
    orgs = find(c==tcg);  %the organisms of current cluster
    x = xg(orgs);  %x coords
    Mx = max(x);  mx = min(x);  %get min and max x coords
    cluster_spans(C,1) = ((Mx-mx)/Lx);
    y = yg(orgs);  %y coords
    My = max(y);  my = min(y);
    cluster_spans(C,2) = ((My-my)/Ly);
  end %c
end %gen
end