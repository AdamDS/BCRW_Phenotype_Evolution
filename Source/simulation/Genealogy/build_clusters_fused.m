%% build_clusters_fused.m
% function [num_clusters_fused,clusters_fused] = build_clusters_fused(base_name,run_name)
% This function is a historical approach to cluster lineages, that is,
% its focus is on "where" clusters came "from". This is like looking at the
% lines of ancestory composing each cluster and therefore is a fusion of
% past clusters.
% Two data vectors are saved, num_clusters_fused & clusters_fused. 
% num_clusters_fused tells the number of parent clusters to each offspring 
% cluster, so it's length is sum(population(2:NGEN)) to index the offspring 
% information. clusters_fused tells which parent clusters go to each offspring
% cluster, so each set of parent clusters is grouped by length of 
% num_clusters_fused(x), in order of offspring cluster id by generation.
% This function uses population, parents, num_clusters, and trace_clusters. 
% It also considers the cluster seed limit value (limit), so it will name the 
% mat files accordingly.
%
% if num_clusters_fused(x) = 0, then impossible
% if num_clusters_fused(x) = 1, then pure cluster from parent's generation
% if num_clusters_fused(x) > 1, then parent clusters mixed (converge)
function [num_clusters_fused,clusters_fused] = build_clusters_fused(...
  base_name,run,dir_name), 
global SIMOPTS;
this_function = 'build_clusters_fused';

new_dir_name = split_cd(dir_name,run,SIMOPTS.split,1,0);
run_name = int2str(run);
clus_name = cluster_name(base_name);

fprintf(' %s for %s \n',this_function,[clus_name run_name]);

num_clusters_fused = [];  clusters_fused = [];

if ~mat_exist([new_dir_name 'clusters_fused_' clus_name run_name]) && ...
   ~mat_exist([new_dir_name 'num_clusters_fused_' clus_name run_name]) || ...
   SIMOPTS.write_over,  
[tc,go,error] = try_catch_load([new_dir_name 'trace_cluster_' clus_name run_name],1,1);
if go, [nc,go,error] = try_catch_load([new_dir_name 'num_clusters_' clus_name run_name],1,1);
if go, [par,go,error] = try_catch_load([new_dir_name 'parents_' base_name run_name],1,1);
if go, [pop,go,error] = try_catch_load([new_dir_name 'population_' base_name run_name],1,1);
if go,             
  trace_cluster = tc.trace_cluster;  clear tc
  num_clusters = nc.num_clusters;  clear nc
  parents = par.parents;  clear par
  population = pop.population;  clear pop

  % begin debug level
  ou = 0; ov = 0; pu = 0; pv = 0; % Generational indices for offspring and parents
  ocu = 0;  ocv = 0;  pcu = 0;  pcv = 0;  % Generational indices for offspring and parent clusters
  ngen = length(find(population));
  IPOP = population(1);
  clusters_fused = [];  %the parent clusters to each offspring cluster
  num_clusters_fused = zeros(sum(num_clusters(2:ngen)),1);  %the number of parent clusters to each offspring cluster
  for gen = 1:(ngen-1), %for each parent generation
    script_gen_update(this_function,gen,base_name,run_name);
    ou = ov +1; ov = sum(population(2:(gen+1)));  %offspring indices
    pu = pv +1; pv = sum(population(1:gen));  %parent indices
    ocu = ocv +1; ocv = sum(num_clusters(2:gen+1)); %offspring cluster indices
    pcu = pcv +1; pcv = sum(num_clusters(1:gen)); %parent cluster indices

    par = parents(ou:ov,:); %parents (from the organisms in gen) of offspring in gen+1
    offspring_clusters = trace_cluster((pv+1):(ov+IPOP)); 
    parent_clusters = trace_cluster(pu:pv);
    for oc = 1:num_clusters(gen+1), %for each offspring cluster in gen+1
      offspring_of_clusters = find(offspring_clusters==oc); %the offspring in cluster oc
      %parents of the offspring in cluster oc
      parents_of_offspring_of_clusters = unique([par(offspring_of_clusters,1); ...
                                                 par(offspring_of_clusters,2)]);
      %clusters the parents were in
      clusters_of_parents = unique(parent_clusters(parents_of_offspring_of_clusters));
      %total number of parent clusters which fused to form offspring cluster oc
      num_clusters_fused(ocu+oc-1) = length(clusters_of_parents);
      %list of the parent clusters which fused for each offspring cluster
      clusters_fused = [clusters_fused; clusters_of_parents];
    end
  end
  % end debug level
  save([new_dir_name 'clusters_fused_' clus_name run_name],'clusters_fused');
  save([new_dir_name 'num_clusters_fused_' clus_name run_name],'num_clusters_fused');
end %population
end %parents
end %num_clusters
end %trace_cluster
end %exists
end %function