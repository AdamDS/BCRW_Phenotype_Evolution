%% Clustering.m
% function [times] = Clustering()
% Determines trace_cluster_seed, seed_distances, num_clusters, trace_clusters, orgsncluster, 
% centroids of clusters and cluster_diversity given the necessary parameter options in SIMOPTS.
% The output gives the time to complete all clustering functions designated.
function [times] = Clustering(),  
global SIMOPTS;
times = zeros(length(SIMOPTS.SIMS),1);
[base_name,dir_name] = NameAndCD(SIMOPTS.make_dir,SIMOPTS.do_cd);
[clus_name] = cluster_name(base_name);
i = 0;  tic;
for run = SIMOPTS.SIMS, 
  if ~SIMOPTS.write_over, fprintf(['Attempting clustering ' clus_name int2str(run) '\n']);  end
  i = i +1;
  SIMOPTS.run = run;
  
  if SIMOPTS.reproduction==1 && SIMOPTS.do_cluster_seeds,
    [~] = build_cluster_seeds(base_name,run,dir_name);  end
  if SIMOPTS.do_build_clusters, 
    [~,~,~] = build_clusters(base_name,run,dir_name); end
  if SIMOPTS.do_locate_clusters,  
    [~,~,~] = locate_clusters(base_name,run,dir_name);  end
  if SIMOPTS.do_build_gyration_radii, 
    [~] = build_gyration_radii(base_name,run,dir_name); end
%     if SIMOPTS.fix, [~] = fix_gyration_radii(base_name,run,dir_name); 
%     else, [~] = build_gyration_radii(base_name,run,dir_name); end, end
  if SIMOPTS.do_correlation_lengths,  
    [~] = build_correlation_lengths(base_name,run,dir_name);  end
  if SIMOPTS.do_build_diameters,  
    if SIMOPTS.pool_size>0, [~] = CPUPAR_build_cluster_diameters(base_name,run_name); 
    elseif SIMOPTS.pool_size<0, [~] = GPUPAR_build_cluster_diameters(base_name,run_name); 
    else, [~] = build_cluster_diameters(base_name,run,dir_name); end
    if SIMOPTS.fix, [~] = fix_cluster_diameters(base_name,run_name);  end
  end
  times(i) = toc;
end
end