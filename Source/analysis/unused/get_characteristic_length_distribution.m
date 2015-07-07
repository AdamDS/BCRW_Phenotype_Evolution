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
CD = [];  CD_finite = [];
tic;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  i = find(mu==mutability);
  for run = SIMS
    run_name = int2str(run);
    exp_name = [base_name run_name];
    if do_characteristic_diameters==1, 
      cd_name = ['cluster_diameters_' exp_name];
      go = 1; [cd,go] = try_catch_load(cd_name,go);
      if go==1
        cluster_diameters = cd.cluster_diameters;  clear cd
        CD = [CD; cluster_diameters];
        CD_finite = [CD_finite; cluster_diameters(cluster_diameters~=Inf)];
      end
    end
  end
end
end
end
toc;

