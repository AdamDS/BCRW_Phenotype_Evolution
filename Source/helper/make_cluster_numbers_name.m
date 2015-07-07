function [cn_name] = make_cluster_numbers_name(dir_name,clus_name,SIMS);
cn_name = [];
try,  
  cn_name = [dir_name 'cluster_numbers_' clus_name(1:end-1) ...
        int2str(SIMS(1)) '_' int2str(SIMS(end) '_' int2str(length(SIMS)))];
catch error,  
  fprintf('%s \n',error);
end
end