% if num_cluster_fusions(x) = 0, then impossible
% if num_cluster_fusions(x) = 1, then pure cluster from parent's generation
% if num_cluster_fusions(x) > 1, then parent clusters mixed (converge)

% if num_clusters_produced(x) = 0, then parent cluster died out
% if num_clusters_produced(x) = 1, then parent cluster did not split
% if num_clusters_produced(x) > 1, then parent clusters split (diverge)
this_script = 'get_cluster_lineage_hists';
fprintf([this_script '\n']);
global SIMOPTS;
legend_colors = ['bkg'];  lc = 0;
% mutability = [0.2 0.38 0.4 0.42 0.48:0.02:0.52 1.20]; 
b = 0;
zm = zeros(length(mutability),1);
lshncfs = zm;
lshncps = zm;
all_shncfs = [];
all_shncps = [];
all_xshncfs = [];
all_xshncps = [];
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  lc = lc +1;
  b = b +1;
  hncfs = [];
  hncps = [];
  for run = SIMS
    run_name = int2str(run);
    load(make_data_name('num_clusters',base_name,run_name,0));
    load(make_data_name('num_clusters_fused',base_name,run_name,0));
    load(make_data_name('num_clusters_produced',base_name,run_name,0));
    hncf = hist(num_clusters_fused,max(num_clusters_fused));
    hncfs = cat_row(hncfs,hncf);
    hncp = hist(num_clusters_produced,max(num_clusters_produced));
    hncps = cat_row(hncps,hncp);
  end
  shncfs = sum(hncfs,1)/length(SIMS);
  shncps = sum(hncps,1)/length(SIMS);
  figure(82000000 +mu*1000 +limit +100000*reproduction), 
  semilogy(shncfs,'b');
  hold on;
  semilogy(0:(length(shncps)-1),shncps,'k');
  title(['Cluster branching ' make_title_name(make_data_name('',base_name,'',0),'')]);
  xlabel('branching');
  ylabel('<occurrences>');
  legend('parent clusters fused into offspring clusters',...
    'offspring clusters produced from parent clusters');

  lshncfs(b) = length(shncfs);
  lshncps(b) = length(shncps);
  all_shncfs = cat_row(all_shncfs,shncfs);
  all_shncps = cat_row(all_shncps,shncps);
  all_xshncfs = cat_row(all_xshncfs,[1:length(shncfs)]);
  all_xshncps = cat_row(all_xshncps,[0:(length(shncps)-1)]);
end
end
end
figure(83000000 +limit +100000*reproduction),
plot(mutability,lshncfs,'bx'); hold on; %length of sum of hist of num clusters fused
plot(mutability,lshncps,'kx'); hold off; %length of sum of hist of num clusters produced
title(make_title_name(generalize_base_name(base_name),''));
xlabel('\mu');  ylabel('greatest interaction');
legend('fused','produced');

figure(84000000 +limit +100000*reproduction), 
mesh(all_shncfs);
yTickLabel(int2str(100*mutability))