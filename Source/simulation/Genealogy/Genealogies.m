%% Genealogies.m
% function [times] = Genealogies()
% Determines (num_descendants), num_descendant_clusters, descendant_clusters,
% num_clusters_fused, clusters_fused, num_clusters_produced, and clusters_produced
% given the necessary parameter options in SIMOPTS.
% The output gives the time to complete each run round (so length is that of SIMS) 
% of genealogy functions.
function [times] = Genealogies(), 
global SIMOPTS;
[base_name,dir_name] = NameAndCD(SIMOPTS.make_dir,SIMOPTS.do_cd);
i = 0;  tic;
for run = SIMOPTS.SIMS, 
  if SIMOPTS.write_over==0, fprintf(['Attempting genealogies ' base_name int2str(run) '\n']); end
  i = i +1;
  SIMOPTS.run = run;
  
  if SIMOPTS.do_indiv_lineage,  
    [~] = build_num_descendants(base_name,run,dir_name);  end
  if SIMOPTS.do_indiv_cluster_lineage,  
    [~,~] = build_descendant_clusters(base_name,run,dir_name); clear ndc dc;  end
  if SIMOPTS.do_cluster_lineage,  
    [~,~] = build_clusters_fused(base_name,run,dir_name); clear ncf cf; 
    [~,~] = build_clusters_produced(base_name,run,dir_name);  clear ncp cp; end
  times(i) = toc;
end
end