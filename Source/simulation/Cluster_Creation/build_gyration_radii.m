%% get_gyration_radii.m
%
% Builds the data vector, cluster_diameters, for a set of data. The
% vector is organized like a centroids data vector, in that the order of
% the elements are blocked off by generation and by cluster identity order.
%
% FIXED ON 22 AUG 2013; ALL gyration_radii DATA BEFORE NEEDS TO BE FIXED!!!
% old: radii(cluster) = sqrt(mean((sum((these_coords-these_cents).^2,2))./orgs));
% new: radii(cluster) = sqrt(mean(sum((these_coords-these_cents).^2,2)));
%
% -ADS 10*4*12
function [gyration_radii] = build_gyration_radii(base_name,run,dir_name), 
global SIMOPTS;
this_function = 'build_gyration_radii';

new_dir_name = split_cd(dir_name,run,SIMOPTS.split,1,0);
run_name = int2str(run);
[clus_name] = cluster_name(base_name);

print_function(this_function,[clus_name run_name]);

gon = 0;  goo = 0;
limit = SIMOPTS.limit;
gyration_radii = [];

if ~mat_exist([new_dir_name 'gyration_radii_' clus_name run_name]) || SIMOPTS.write_over, 
[p,go,error] = exist_load([new_dir_name 'population_' base_name run_name],1,1);
if go, [tc,go,error] = exist_load([new_dir_name 'trace_cluster_' clus_name run_name],1,1);
if go, [nc,go,error] = exist_load([new_dir_name 'num_clusters_' clus_name run_name],1,1);
if go, [tx,go,error] = exist_load([new_dir_name 'trace_x_' base_name run_name],1,1);
if go, [ty,go,error] = exist_load([new_dir_name 'trace_y_' base_name run_name],1,1);
if go, [cx,go,error] = exist_load([new_dir_name 'centroid_x_' clus_name run_name],1,1);
if go, [cy,go,error] = exist_load([new_dir_name 'centroid_y_' clus_name run_name],1,1);
if go, [org,go,error] = exist_load([new_dir_name 'orgsnclusters_' clus_name run_name],1,1);
if go,  
  fprintf([clus_name run_name '\n']);
  trace_cluster = tc.trace_cluster;  clear tc
  num_clusters = nc.num_clusters;  clear nc
  trace_x = tx.trace_x; clear tx
  trace_y = ty.trace_y; clear ty
  population = p.population;  clear pop
  centroid_x = cx.centroid_x; clear cx
  centroid_y = cy.centroid_y; clear cy
  orgsnclusters = org.orgsnclusters;  clear org

  %initialize characteristic length holders
  gyration_radii = zeros(sum(num_clusters),1); %indexed by j
  ngen = length(population(population>=limit));
  for gen = 1:ngen, 
    script_gen_update(this_function,gen,base_name,run_name);
    v = sum(population(1:gen)); u = v -population(gen) +1;
    cv = sum(num_clusters(1:gen));  cu = cv -num_clusters(gen) +1;
    %get the cluster assignments of organisms in this generation
    tc_of_gen = trace_cluster(u:v);
    %get the centroids of each cluster
    cents_of_gen = [centroid_x(cu:cv) centroid_y(cu:cv)];
    %get the coordinates of organisms in this generation
    coords_of_gen = [trace_x(u:v) trace_y(u:v)];
    orgs_of_gen = orgsnclusters(cu:cv);
    radii = zeros(num_clusters(gen),1);
    %for each cluster
    for cluster = 1:num_clusters(gen),  
      %get organisms of the same cluster
      these_orgs = find(tc_of_gen==cluster);
      N = length(these_orgs);
      %get x & y coords of only those in the same cluster
      these_coords = coords_of_gen(these_orgs,:);
      these_cents = ones(size(these_coords,1),1)*cents_of_gen(cluster,:);
      orgs = orgs_of_gen(cluster);
      %gyration radius for distances between center of cluster mass and
      %each organism (Lesne p15 assuming integral A(xbar) = orgs)
%       radii2(cluster) = sqrt(mean((sum((these_coords-these_cents).^2,2))./orgs));
%FIXED ON 22 AUG 2013; ALL gyration_radii DATA BEFORE NEEDS TO BE FIXED!!!
      radii(cluster) = sqrt(mean(sum((these_coords-these_cents).^2,2)));
%% start of debugging
%%
    end
    gyration_radii(cu:cv) = radii;
%     gyration_radii2(cu:cv) = radii2;
  end
  save([new_dir_name 'gyration_radii_' clus_name run_name],'gyration_radii');
%   save([new_dir_name 'gyration_radii2_' clus_name run_name],'gyration_radii2');
%% 
%% last part of debugging
end %orgsnclusters
end %centroid_y
end %centroid_x
end %trace_y
end %trace_x
end %num_clusters
end %trace_cluster
end %population

end %exists
%   end %for SIMS
% end %for mutability
% end %for death_max
% end %for overpop
end %function