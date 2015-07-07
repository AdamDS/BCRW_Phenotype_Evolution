% if num_cluster_fusions(x) = 0, then impossible
% if num_cluster_fusions(x) = 1, then pure cluster from parent's generation
% if num_cluster_fusions(x) > 1, then parent clusters mixed (converge)

% if num_clusters_produced(x) = 0, then parent cluster died out
% if num_clusters_produced(x) = 1, then parent cluster did not split
% if num_clusters_produced(x) > 1, then parent clusters split (diverge)
this_script = 'get_lineage_stats';
fprintf([this_script '\n']);
global SIMOPTS;
legend_colors = ['bkg'];  lc = 0;
% mutability = [0.2 0.38 0.4 0.42 0.48:0.02:0.52 1.20]; 
b = 0;
zm = zeros(length(mutability),1);
lshncfs = zm;
lshncps = zm;
extinctions = zm;
std_e = zm;
pures = zm;
no_splits = zm;
convergences = zm;
divergences = zm;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  lc = lc +1;
  b = b +1;
  hncfs = [];
  hncps = [];
  extinction = [];
  pure = [];
  no_split = [];
  convergence = [];
  divergence = [];
  for run = SIMS
    run_name = int2str(run);
    load(make_data_name('num_clusters',base_name,run_name,0));
    load(make_data_name('num_clusters_fused',base_name,run_name,0));
    load(make_data_name('num_clusters_produced',base_name,run_name,0));
    hncf = hist(num_clusters_fused,max(num_clusters_fused));
    hncfs = cat_row(hncfs,hncf);
    hncp = hist(num_clusters_produced,max(num_clusters_produced));
    hncps = cat_row(hncps,hncp);
    extinction = [extinction; hncp(1)];
    pure = [pure; hncf(1)];
    no_split = [no_split; hncp(2)];
    convergence = [convergence; hncf(2:end)'];
    divergence = [divergence; hncp(3:end)'];
  end
  extinctions(b) = mean(extinction);
  std_e(b) = std(extinction);
  pures(b) = mean(pure);
  std_p(b) = std(pure);
  no_splits(b) = mean(no_split);
  std_ns(b) = std(no_split);
  convergences(b) = mean(convergence);
  std_c(b) = std(convergence);
  divergences(b) = mean(divergence);
  std_d(b) = std(divergence);
  shncfs = sum(hncfs,1)/length(SIMS);
  shncps = sum(hncps,1)/length(SIMS);
  lshncfs(b) = length(shncfs);
  lshncps(b) = length(shncps);
end
end
end
if do_cl_stat_plots==1
  figure(91000000 +100000*reproduction +1000*mu +limit), errorbar(mutability,extinctions,std_e);
  xlim([min(mutability)-dmu max(mutability)+dmu]);  ylim([0 max(extinctions)+max(std_e)]);
  title(make_title_name(generalize_base_name(base_name),''));
  xlabel('\mu');  ylabel('extinctions');
  
  figure(92000000 +100000*reproduction +1000*mu +limit), errorbar(mutability,pures,std_p);
  xlim([min(mutability)-dmu max(mutability)+dmu]);  ylim([0 max(pures)+max(std_p)]);
  title(make_title_name(generalize_base_name(base_name),''));
  xlabel('\mu');  ylabel('pures');
  
  figure(93000000 +100000*reproduction +1000*mu +limit), errorbar(mutability,no_splits,std_ns);
  xlim([min(mutability)-dmu max(mutability)+dmu]);  ylim([0 max(no_splits)+max(std_ns)]);
  title(make_title_name(generalize_base_name(base_name),''));
  xlabel('\mu');  ylabel('no\_splits');
  
  figure(94000000 +100000*reproduction +1000*mu +limit), errorbar(mutability,convergences,std_c);
  xlim([min(mutability)-dmu max(mutability)+dmu]);  ylim([0 max(convergences)+max(std_c)]);
  title(make_title_name(generalize_base_name(base_name),''));
  xlabel('\mu');  ylabel('convergences');
  
  figure(95000000 +100000*reproduction +1000*mu +limit), errorbar(mutability,divergences,std_d);
  xlim([min(mutability)-dmu max(mutability)+dmu]);  ylim([0 max(divergences)+max(std_d)]);
  title(make_title_name(generalize_base_name(base_name),''));
  xlabel('\mu');  ylabel('divergences');

  figure(96000000 +100000*reproduction +1000*mu +limit), 
  EXS = extinctions./(extinctions +pures +no_splits +convergences +divergences);
  PRS = pures./(extinctions +pures +no_splits +convergences +divergences);
  NSS = no_splits./(extinctions +pures +no_splits +convergences +divergences);
  CONS = convergences./(extinctions +pures +no_splits +convergences +divergences);
  DIVS = divergences./(extinctions +pures +no_splits +convergences +divergences);
  plot(mutability,EXS,'kx',mutability,PRS,'bd',mutability,NSS,'ks',mutability,CONS,'b+',mutability,DIVS,'k*');
  xlim([min(mutability)-dmu max(mutability)+dmu]);  ylim([0 1]);
  title(make_title_name(generalize_base_name(base_name),''));
  xlabel('\mu');  ylabel('proportion');
  legend('extinctions','pures','no\_splits','convergences','divergences');
end