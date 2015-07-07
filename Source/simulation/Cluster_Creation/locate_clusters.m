%% locate_clusters.m
% function [centroid_x,centroid_y,cluster_diversity] = locate_clusters(base_name,run,dir_name), 
%
% -ADS modified from ND
%
% function [centroid_x,centroid_y,cluster_diversity,radii_gyration] = locate_clusters(...
%   base_name,run,dir_name), 
%
function [centroid_x,centroid_y,cluster_diversity] = locate_clusters(base_name,run,dir_name), 
global SIMOPTS;
this_function = 'locate_clusters';

new_dir_name = split_cd(dir_name,run,SIMOPTS.split,1,0);
run_name = int2str(run);

centroid_x = [];  centroid_y = [];  cluster_diversity = [];
[clus_name] = cluster_name(base_name);

print_function(this_function,[clus_name run_name]);

if ~mat_exist([new_dir_name 'centroid_x_' clus_name run_name]) || ...
   ~mat_exist([new_dir_name 'centroid_y_' clus_name run_name]) || ...
   ~mat_exist([new_dir_name 'cluster_diversity_' clus_name run_name]) || ...
   SIMOPTS.write_over, 

  [p,go] = try_catch_load([new_dir_name 'population_' base_name run_name],1,1);  
  if go,  [tx,go] = try_catch_load([new_dir_name 'trace_x_' base_name run_name],1,1);
  if go,  [ty,go] = try_catch_load([new_dir_name 'trace_y_' base_name run_name],1,1);
  if go,  [tc,go] = try_catch_load([new_dir_name 'trace_cluster_' clus_name run_name],1,1);
  if go,  [nc,go] = try_catch_load([new_dir_name 'num_clusters_' clus_name run_name],1,1);
  if go,  

    population = p.population;  clear p
    trace_x = tx.trace_x; clear tx
    trace_y = ty.trace_y; clear ty
    trace_cluster = tc.trace_cluster; clear tc
    num_clusters = nc.num_clusters; clear nc

    NGEN = length(find(population));
    centroid_x = zeros(sum(num_clusters),1);
    centroid_y = zeros(sum(num_clusters),1);
    cluster_diversity = zeros(sum(num_clusters),1);
    radii_gyration = zeros(sum(num_clusters),1);
    
    low_p = 0; high_p = 0; low_nc = 0; high_nc = 0;

    for gen = 1:NGEN, 
      low_p = high_p +1; high_p = sum(population(1:gen)); 
      indiv_x = trace_x(low_p:high_p);  indiv_y = trace_y(low_p:high_p);
      low_nc = high_nc +1; high_nc = sum(num_clusters(1:gen));
      cluster_labels = trace_cluster(low_p:high_p); 
      for this_cluster = 1:num_clusters(gen), 
        these_indivs = find(cluster_labels==this_cluster);  %indivs within cluster
        centroid_x(low_nc+this_cluster-1) = mean(indiv_x(these_indivs));  %calculate centroid_x
        centroid_y(low_nc+this_cluster-1) = mean(indiv_y(these_indivs));  %calculate centroid_y
        %calculate cluster_diversity
        drs = pdist([indiv_x(these_indivs),indiv_y(these_indivs)]);
        cluster_diversity(low_nc+this_cluster-1) = mean(drs);
        d = 0; k = 0;
        for i=1:length(these_indivs), %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          for j=1:length(these_indivs), %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if i~=j
              k = k +1;
              dx = [indiv_x(these_indivs(i))-indiv_x(these_indivs(j))].^2;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              dy = [indiv_y(these_indivs(i))-indiv_y(these_indivs(j))].^2;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              d(k) = dx+dy;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
          end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        radii_gyration(low_nc+this_cluster-1) = sqrt(sum(d)/2)/length(these_indivs);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      end
    end

    save([new_dir_name 'centroid_x_' clus_name run_name],'centroid_x');
    save([new_dir_name 'centroid_y_' clus_name run_name],'centroid_y');
    save([new_dir_name 'cluster_diversity_' clus_name run_name],'cluster_diversity');
    save([new_dir_name 'radii_gyration_' clus_name run_name],'radii_gyration');
  end %num_clusters
  end %trace_cluster
  end %trace_y
  end %trace_x
  end %population
end %exists
end %function

%% Old Code
% cstd_x = zeros(sum(num_clusters),1);
% cskew_x = zeros(sum(num_clusters),1);
% cstd_y = zeros(sum(num_clusters),1);
% cskew_y = zeros(sum(num_clusters),1);
%     cstd_x(low_nc+this_cluster-1) = std(indiv_x(these_indivs));
%     cskew_x(low_nc+this_cluster-1) = skewness(indiv_x(these_indivs));
%     cstd_y(low_nc+this_cluster-1) = std(indiv_y(these_indivs));
%     cskew_y(low_nc+this_cluster-1) = skewness(indiv_y(these_indivs));