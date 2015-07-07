%% movefile_data_set.m
% Inputs:
% old_name = string of current base name of data set
% new_name = string of desired base name of data set
% run_name = string of run integer
%
%%
function [] = movefile_data_set(old_name,new_name,run_name,opt),  
if opt, 
  if exist(['population_' old_name run_name '.mat'])==2
    movefile(['population_' old_name run_name '.mat'],['population_' new_name run_name '.mat']);
  end

  if exist(['trace_x_' old_name run_name '.mat'])==2
    movefile(['trace_x_' old_name run_name '.mat'],['trace_x_' new_name run_name '.mat']);
  end

  if exist(['trace_y_' old_name run_name '.mat'])==2
    movefile(['trace_y_' old_name run_name '.mat'],['trace_y_' new_name run_name '.mat']);
  end

  if exist(['kills_' old_name run_name '.mat'])==2
    movefile(['kills_' old_name run_name '.mat'],['kills_' new_name run_name '.mat']);
  end

  if exist(['rivalries_' old_name run_name '.mat'])==2
    movefile(['rivalries_' old_name run_name '.mat'],['rivalries_' new_name run_name '.mat']);
  end

  if exist(['parents_' old_name run_name '.mat'])==2
    movefile(['parents_' old_name run_name '.mat'],['parents_' new_name run_name '.mat']);
  end

  if exist(['trace_noise_' old_name run_name '.mat'])==2
    movefile(['trace_noise_' old_name run_name '.mat'],['trace_noise_' new_name run_name '.mat']);
  end

  if exist(['trace_cluster_' old_name run_name '.mat'])==2
    movefile(['trace_cluster_' old_name run_name '.mat'],['trace_cluster_' new_name run_name '.mat']);
  end

  if exist(['trace_cluster_seed_' old_name run_name '.mat'])==2
    movefile(['trace_cluster_seed_' old_name run_name '.mat'],['trace_cluster_seed_' new_name run_name '.mat']);
  end

  if exist(['seed_distances_' old_name run_name '.mat'])==2
    movefile(['seed_distances_' old_name run_name '.mat'],['seed_distances_' new_name run_name '.mat']);
  end

  if exist(['num_clusters_' old_name run_name '.mat'])==2
    movefile(['num_clusters_' old_name run_name '.mat'],['num_clusters_' new_name run_name '.mat']);
  end

  if exist(['centroidx_' old_name run_name '.mat'])==2
    load(['centroidx_' old_name run_name '.mat']);
    centroid_x = centroidx; clear centroidx
    save(['centroid_x_' new_name run_name '.mat'],'centroid_x');
  elseif exist(['centroid_x_' old_name run_name '.mat'])==2
    movefile(['centroid_x_' old_name run_name '.mat'],['centroid_x_' new_name run_name '.mat']);
  end

  if exist(['centroidy_' old_name run_name '.mat'])==2
    load(['centroidy_' old_name run_name '.mat']);
    centroid_y = centroidy; clear centroidy
    save(['centroid_y_' new_name run_name '.mat'],'centroid_y');
  elseif exist(['centroid_y_' old_name run_name '.mat'])==2
    movefile(['centroid_y_' old_name run_name '.mat'],['centroid_y_' new_name run_name '.mat']);
  end

  if exist(['dens_clusters_' old_name run_name '.mat'])==2
    load(['dens_clusters_' old_name run_name '.mat']);
    cluster_diversity = dens_clusters;  clear dens_clusters
    save(['cluster_diversity_' new_name run_name '.mat'],'cluster_diversity');
  elseif exist(['cluster_diversity_' old_name run_name '.mat'])==2
    movefile(['cluster_diversity_' old_name run_name '.mat'],['cluster_diversity_' new_name run_name '.mat']);
  end

  if exist(['orgsnclusters_' old_name run_name '.mat'])==2
    movefile(['orgsnclusters_' old_name run_name '.mat'],['orgsnclusters_' new_name run_name '.mat']);
  end

  if exist(['cluster_coefficient_' old_name run_name '.mat'])==2
    movefile(['cluster_coefficient_' old_name run_name '.mat'],['cluster_coefficient_' new_name run_name '.mat']);
  end

  if exist(['in_degree_' old_name run_name '.mat'])==2
    movefile(['in_degree_' old_name run_name '.mat'],['in_degree_' new_name run_name '.mat']);
  end

  if exist(['fused_clusters_' old_name run_name '.mat'])==2
    movefile(['fused_clusters_' old_name run_name '.mat'],['fused_clusters_' new_name run_name '.mat']);
  end
  
else, %change by new directory
  old_slash = find(old_name=='\');
  o1 = old_name(1:old_slash(end));  %G:...\ (to last slash of directory)
  o2 = old_name(old_slash(end)+1:end);  %Bacterial...2_Flatscape_4 (from last slash of directory)
  new_slash = find(new_name=='\');
  n1 = new_name(1:new_slash(end));  %G:...\ (to last slash of directory)
  n2 = new_name(new_slash(end)+1:end);  %Bacterial...2_Flatscape_4 (from last slash of directory)

  if exist([o1 'population_' o2 run_name '.mat'])==2
    movefile([o1 'population_' o2 run_name '.mat'],[n1 'population_' n2 run_name '.mat']);
  end

  if exist([o1 'trace_x_' o2 run_name '.mat'])==2
    movefile([o1 'trace_x_' o2 run_name '.mat'],[n1 'trace_x_' n2 run_name '.mat']);
  end

  if exist([o1 'trace_y_' o2 run_name '.mat'])==2
    movefile([o1 'trace_y_' o2 run_name '.mat'],[n1 'trace_y_' n2 run_name '.mat']);
  end

  if exist([o1 'kills_' o2 run_name '.mat'])==2
    movefile([o1 'kills_' o2 run_name '.mat'],[n1 'kills_' n2 run_name '.mat']);
  end

  if exist([o1 'rivalries_' o2 run_name '.mat'])==2
    movefile([o1 'rivalries_' o2 run_name '.mat'],[n1 'rivalries_' n2 run_name '.mat']);
  end

  if exist([o1 'parents_' o2 run_name '.mat'])==2
    movefile([o1 'parents_' o2 run_name '.mat'],[n1 'parents_' n2 run_name '.mat']);
  end

  if exist([o1 'trace_noise_' o2 run_name '.mat'])==2
    movefile([o1 'trace_noise_' o2 run_name '.mat'],[n1 'trace_noise_' n2 run_name '.mat']);
  end

  if exist([o1 'trace_cluster_' o2 run_name '.mat'])==2
    movefile([o1 'trace_cluster_' o2 run_name '.mat'],[n1 'trace_cluster_' n2 run_name '.mat']);
  end

  if exist([o1 'trace_cluster_seed_' o2 run_name '.mat'])==2
    movefile([o1 'trace_cluster_seed_' o2 run_name '.mat'],[n1 'trace_cluster_seed_' n2 run_name '.mat']);
  end

  if exist([o1 'seed_distances_' o2 run_name '.mat'])==2
    movefile([o1 'seed_distances_' o2 run_name '.mat'],[n1 'seed_distances_' n2 run_name '.mat']);
  end

  if exist([o1 'num_clusters_' o2 run_name '.mat'])==2
    movefile([o1 'num_clusters_' o2 run_name '.mat'],[n1 'num_clusters_' n2 run_name '.mat']);
  end

  if exist([o1 'centroidx_' o2 run_name '.mat'])==2
    load([o1 'centroidx_' o2 run_name '.mat']);
    centroid_x = centroidx; clear centroidx
    save([o1 'centroid_x_' n2 run_name '.mat'],'centroid_x');
  elseif exist([o1 'centroid_x_' o2 run_name '.mat'])==2
    movefile([o1 'centroid_x_' o2 run_name '.mat'],[n1 'centroid_x_' n2 run_name '.mat']);
  end

  if exist([o1 'centroidy_' o2 run_name '.mat'])==2
    load([o1 'centroidy_' o2 run_name '.mat']);
    centroid_y = centroidy; clear centroidy
    save([o1 'centroid_y_' n2 run_name '.mat'],'centroid_y');
  elseif exist([o1 'centroid_y_' o2 run_name '.mat'])==2
    movefile([o1 'centroid_y_' o2 run_name '.mat'],[n1 'centroid_y_' n2 run_name '.mat']);
  end

  if exist([o1 'dens_clusters_' o2 run_name '.mat'])==2
    load([o1 'dens_clusters_' o2 run_name '.mat']);
    cluster_diversity = dens_clusters;  clear dens_clusters
    save([o1 'cluster_diversity_' n2 run_name '.mat'],'cluster_diversity');
  elseif exist([o1 'cluster_diversity_' o2 run_name '.mat'])==2
    movefile([o1 'cluster_diversity_' o2 run_name '.mat'],[n1 'cluster_diversity_' n2 run_name '.mat']);
  end

  if exist([o1 'orgsnclusters_' o2 run_name '.mat'])==2
    movefile([o1 'orgsnclusters_' o2 run_name '.mat'],[n1 'orgsnclusters_' n2 run_name '.mat']);
  end

  if exist([o1 'cluster_coefficient_' o2 run_name '.mat'])==2
    movefile([o1 'cluster_coefficient_' o2 run_name '.mat'],[n1 'cluster_coefficient_' n2 run_name '.mat']);
  end

  if exist([o1 'in_degree_' o2 run_name '.mat'])==2
    movefile([o1 'in_degree_' o2 run_name '.mat'],[n1 'in_degree_' n2 run_name '.mat']);
  end

  if exist([o1 'fused_clusters_' o2 run_name '.mat'])==2
    movefile([o1 'fused_clusters_' o2 run_name '.mat'],[n1 'fused_clusters_' n2 run_name '.mat']);
  end
end