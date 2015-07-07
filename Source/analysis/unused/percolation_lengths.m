%percolation_lengths.m
% ????need to define percolation_limit (0 to 1, represents percentage of
% landscape length largest cluster occupies)
global SIMOPTS;
avg_max_perc_x = zeros(length(mutability),length(SIMS));
avg_max_perc_y = zeros(length(mutability),length(SIMS));
AVG_MAX_PERC_X = zeros(length(mutability),1);
STD_MAX_PERC_X = zeros(length(mutability),1);
AVG_MAX_PERC_Y = zeros(length(mutability),1);
STD_MAX_PERC_Y = zeros(length(mutability),1);
i = 0;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  i = i +1;
  Lx = ((2*((2*basic_map_size(1)) -1)) -1);
  Ly = ((2*((2*basic_map_size(2)) -1)) -1);
  for run = SIMS
    run_name = int2str(run);
    %name the files to load
    p_name = make_data_name('population',base_name,run_name,0);
    tx_name = make_data_name('trace_x',base_name,run_name,0);
    ty_name = make_data_name('trace_y',base_name,run_name,0);
    nc_name = make_data_name('num_clusters',base_name,run_name,0);
    tc_name = make_data_name('trace_cluster',base_name,run_name,0);
    %load necessary files
    if exist([p_name '.mat'])==2, load(p_name); end
    if exist([tx_name '.mat'])==2,  load(tx_name);  end
    if exist([ty_name '.mat'])==2,  load(ty_name);  end
    if exist([nc_name '.mat'])==2,  load(nc_name);  end
    if exist([tc_name '.mat'])==2,  load(tc_name);  end
    %determine run specific constants
    gens = [(length(num_clusters) -floor(0.1*length(num_clusters))):length(num_clusters)];
    %initialize run data information variables
    max_perc_x = zeros(length(gens),1);
    max_perc_y = zeros(length(gens),1);
    %initialize generation counters
    g = 0;  %g counts the absolute number of generations of interest 
    u = 0;  v = 0;  %u and v determine specific generational data locations
    %loop through actual generations of interest
    for gen = 1:length(population)%gens %gen counts the actual generations of interest
      g = g +1;
      u = v +1; v = sum(population(1:gen));
      perc_x = zeros(num_clusters(gen),1);
      perc_y = zeros(num_clusters(gen),1);
      %loop through clusters of generation gen
      for c = 1:num_clusters(gen)
        orgsofcluster = find(trace_cluster(u:v)==c);
        perc_x(c) = (max(trace_x(orgsofcluster)) -min(trace_x(orgsofcluster))) /Lx;
        perc_y(c) = (max(trace_y(orgsofcluster)) -min(trace_y(orgsofcluster))) /Ly;
      end
      max_perc_x(g) = max(perc_x);
      max_perc_y(g) = max(perc_y);
    end
    pops(1:length(population),run) = population';
    avg_max_perc_x(i,run) = mean(max_perc_x);
    avg_max_perc_y(i,run) = mean(max_perc_y);
  end
  AVG_MAX_PERC_X(i) = mean(avg_max_perc_x(i,:));
  STD_MAX_PERC_X(i) = std(avg_max_perc_x(i,:));
  AVG_MAX_PERC_Y(i) = mean(avg_max_perc_y(i,:));
  STD_MAX_PERC_Y(i) = std(avg_max_perc_y(i,:));
end
end
end

if do_percolation_plot==1
  figure(724000 +reproduction*10 +basic_map_size(1));
  tn = make_title_name(generalize_base_name(base_name),'');
  errorbar(mutability,AVG_MAX_PERC_X,STD_MAX_PERC_X,'xk'); hold on;  
  errorbar(mutability,AVG_MAX_PERC_Y,STD_MAX_PERC_Y,'xb');
  title(tn,'FontSize',16);  xlabel('\mu','FontSize',14);  
  ylabel('<max cluster percolation length> / landscape size','FontSize',14);
%   set(gca,'XTick',[min(mutability):((max(mutability)-min(mutability))/5):max(mutability)],...
%       'Box','on','FontSize',14,'TickDir','in','TickLength',[0.015,0.025]);
%   for p = min(mutability):((max(mutability)-min(mutability))/5):max(mutability),  set(gca,'XTickLabel',{num2str(p)}); end
end