% function [abundances,num_bins] = get_abundances_curve(do_abundances_plot,log_abundance)
global SIMOPTS;
i = 0;  flag = 1;
AVG_AB = [];  avg_abs = [];
num_rel_abundances = zeros(length(mutability));
RELATIVE_ABUNDANCES = [];
RANKS = [];
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  % single simulation analysis
  rel_abundances = [];
  for run = SIMS
    run_name = int2str(run);
    i = i +1;
    load(['num_clusters_' base_name run_name]);
    load(['orgsnclusters_' base_name run_name]);
    ngen = length(find(num_clusters));
    num_rel_abundances(i) = sum(num_clusters);
    [rel_abundance] = sort((orgsnclusters/sum(orgsnclusters)),2,'descend');
    rla = length(rel_abundance);
    rlas = size((rel_abundances),2);
    if rla>rlas
      qwer = rel_abundances;  clear rel_abundances
      rel_abundances = [[qwer zeros(size(qwer,1),rla-rlas)]; rel_abundance];
    else
      rel_abundances = [[rel_abundances]; [rel_abundance zeros(1,rlas-rla)]];
    end
    
  end %end run
  avg_rab_bn = mean(rel_abundances,1); %average over all runs
  std_rab_bn = std(rel_abundances,1);  %std over all runs
  RELATIVE_ABUNDANCES = [RELATIVE_ABUNDANCES [avg_rab_bn; std_rab_bn]]; %#ok<AGROW>
end %end mu
end
end
if record_rank_abundance==1
  save_rank_abundance;
end