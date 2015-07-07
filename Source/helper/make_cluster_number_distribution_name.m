function [cnd_name] = make_cluster_number_distribution_name(dir_name,clus_name,SIMS);
cnd_name = [dir_name 'cluster_number_distribution_' clus_name(1:end-1) ...
      int2str(SIMS(1)) '_' int2str(SIMS(end) '_' int2str(length(SIMS)))];
end