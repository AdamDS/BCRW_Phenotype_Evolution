global SIMOPTS;
fit_level = landscape_heights(1);
i = 0;
avg_nc = zeros(length(mutability),length(SIMS));
min_nc = zeros(length(mutability),length(SIMS));
AVG_NCS = zeros(length(mutability),1);
STD_NCS = zeros(length(mutability),1);
MIN_NCS = zeros(length(mutability),1);
NGENS = zeros(length(mutability),length(SIMS));

AVG_GEN_NCS = zeros(length(mutability),NGEN);
STD_GEN_NCS = zeros(length(mutability),NGEN);

for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir,do_cd);
  gen_ncs = zeros(length(SIMS),NGEN);
  i = i +1;
  for run = SIMS
    new_dir_name = split_cd(dir_name,run,split,make_dir,do_cd);
    run_name = int2str(run);
    [nc,go] = try_catch_load([new_dir_name 'num_clusters_' base_name run_name],1,1);
    if go, 
      num_clusters = nc.num_clusters; clear nc
      g = find(num_clusters);
      if ~by_generation,  
        avg_nc(i,run) = mean(num_clusters);
        min_nc(i,run) = min(num_clusters);
        NGENS(i,run) = length(num_clusters);
      elseif by_generation && ~by_window, 
        gen_ncs(run,1:length(g)) = num_clusters(g);
      elseif by_generation && by_window,  
        j = 0;  overlap = 5;  window = 10;  num_windows = ceil(ngen/overlap); k = -overlap +1;
%         for win = 1:num_windows
%           k = k +overlap;
%           win_ncs(j,run) = num_clusters(k
%         end
      end
    end
  end
  if by_generation, 
    AVG_GEN_NCS(i,:) = mean(gen_ncs,1);
    STD_GEN_NCS(i,:) = std(gen_ncs);
  end
  AVG_NCS(i) = mean(avg_nc(i,:));
  STD_NCS(i) = std(avg_nc(i,:));
  MIN_NCS(i) = min(min_nc(i,:));
end
end
end

if do_num_clusters_plot && ~by_generation,  
  fn18 = 18;
  tn = make_title_name(make_data_name('num_clusters',base_name,'',1),'');
  nz = find(AVG_NCS);
  %errorbar(mutability(nz),AVG_NCS(nz),STD_NCS(nz));
  plot_order_v_control(mutability,AVG_NCS,STD_NCS,fn18,'\mu','num\_clusters','k',base_name);
  title(tn,'FontSize',16);  %xlabel('\mu','FontSize',14);  ylabel('num\_clusters','FontSize',14);
  xlim([min(mutability) max(mutability)]);  %ylim([0 350]); 
%   set(gca,'XTick',[min(bn):((max(bn)-min(bn))/5):max(bn)],...
%       'Box','on','FontSize',14,'TickDir','in','TickLength',[0.015,0.025]);
%   for p = min(bn):((max(bn)-min(bn))/5):max(bn),  set(gca,'XTickLabel',{num2str(p)}); end
end