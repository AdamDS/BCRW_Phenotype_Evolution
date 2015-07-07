%% get_popuations.m *******************************************************
% [avg_pop,min_pop,AVG_POPS,STD_POPS,MIN_POPS] = get_populations(SIMS,...
%   death_max,NGEN,range,IPOP,landscape_movement,landscape_heights,...
%   basic_map_size,flat,shock,rndm,exp_type,reproduct,limit,bn,bi,loaded,load_name,...
%   do_plot)

% function [avg_pop,min_pop,AVG_POPS,STD_POPS,MIN_POPS] = get_populations(SIMS,...
%   death_max,NGEN,range,IPOP,landscape_movement,landscape_heights,...
%   basic_map_size,flat,shock,rndm,exp_type,reproduct,limit,bn,bi,loaded,load_name,...
%   do_plot)
global SIMOPTS;
i = 1;
pncc = zeros(length(mutability),length(SIMS));
tic;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  for run = SIMS
    run_name = int2str(run);
    exp_name = [base_name run_name];
    pop_name = ['population_' exp_name];
    if limit==3
      nc_name = ['num_clusters_' exp_name];
    else
      nc_name = ['num_clusters_' base_name int2str(limit) '_limit_' run_name];
    end
    if exist([pop_name '.mat'])==2 && exist([nc_name '.mat'])==2
      load(pop_name);
      load(nc_name);
      g = find(population>=limit);
      ccpnc = corrcoef(population(g),num_clusters(g));
      pncc(i,run) = ccpnc(1,2);
    end
  end
  i = i +1;
end
end
end
toc;

AVG_PNCC = mean(pncc,2);
for n = 1:length(mutability)
  STD_PNCC(n) = std(pncc(n,:));
end

if do_corr_plot==1
  figure(6549+basic_map_size(1));
  if limit==3
    tn = make_title_name(generalize_base_name(base_name),'');
  else
    tn = [make_title_name(generalize_base_name(base_name),'') '\_' int2str(limit) '\_limit'];
  end

  errorbar(mutability,AVG_PNCC,STD_PNCC,'x');  
  title(['population & num\_clusters correlation ' tn],'FontSize',16);
  xlabel('\mu','FontSize',14);  ylabel('correlation coefficient','FontSize',14);
end