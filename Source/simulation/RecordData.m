%% RecordData.m
%function [] = RecordData(base_name,GENVARS,GENOUTS)
%inputs - base_name (simulation descriptor name), GENVARS (generation variables), &
%GENOUTS (generation data output)
%outputs - none
%
%RecordData is an update of Record4 to account for globalization of SIMOPTS and handling of the
%generation variables & parameters.
%RecordData is set up to work with only the cluster seed information and does nothing about actual
%cluster data. One must run Clustering.m from main_babies.m to generate the cluster data from
%trace_cluster_seed.
%One must also run Genealogies.m from main_babies.m to generate the
%genealogy information.
%
% GENVARS = struct('babies',[],'basic_map',[],'land',[],'run',[],...
%   'exp_name',[],'dir_name',[]);
%
% GENOUTS = struct('population',population,'trace_noise',trace_noise,'trace_x',trace_x,...
%   'trace_y',trace_y,'trace_cluster_seed',trace_cluster_seed,'seed_distances',seed_distances,...
%   'parents',parents,'kills',kills,'rivalries',rivalries,'land',land,'shifted',shifted,...
%   'finished',finished);
%
% UPDATE 14 Aug 2013 - ADS
% Splits data into chunks of SIMOPTS.split by mod(run,split) folders
% No longer changes directories, and instead uses the dir_name from NameAndCD
%
%%
function [] = RecordData(base_name,GENVARS,GENOUTS,dir_name)
global SIMOPTS;

if ~SIMOPTS.only_lt && ~SIMOPTS.only_relax, 
  new_dir_name = split_cd(dir_name,GENVARS.run,SIMOPTS.split,1,0);
end

run_name = int2str(GENVARS.run);
clus_name = cluster_name(base_name);

if ~SIMOPTS.append,
  population = GENOUTS.population;
  gens = find(population>=SIMOPTS.limit);
  len_p = length(gens);
  population = population(gens);
  pop_name = [new_dir_name 'population_' base_name run_name];
  save(pop_name,'population');  clear population
  if SIMOPTS.record,  
    if SIMOPTS.exp_type~=0, 
      trace_noise = GENOUTS.trace_noise;
      tn_name = [new_dir_name 'trace_noise_' base_name run_name];
      save(tn_name,'trace_noise');  clear trace_noise
    end
    trace_x = GENOUTS.trace_x;
    tx_name = [new_dir_name 'trace_x_' base_name run_name];
    save(tx_name,'trace_x');  clear trace_x
    trace_y = GENOUTS.trace_y;
    ty_name = [new_dir_name 'trace_y_' base_name run_name];
    save(ty_name,'trace_y');  clear trace_y
  end
  if SIMOPTS.save_seeds,  
    if SIMOPTS.reproduction~=1, 
      trace_cluster_seed = GENOUTS.trace_cluster_seed;
      tcs_name = [new_dir_name 'trace_cluster_seed_' clus_name run_name];
      save(tcs_name,'trace_cluster_seed');  clear trace_cluster_seed
      seed_distances = GENOUTS.seed_distances;
      sd_name = [new_dir_name 'seed_distances_' clus_name run_name];
      save(sd_name,'seed_distances'); clear seed_distances
    end
  end
  if size(unique(GENOUTS.land),1)~=1, 
    land = GENOUTS.land;
    l_name = [new_dir_name 'land_' base_name run_name];
    save(l_name,'land');  clear land
  end
  if size(GENOUTS.shifted,1)~=0, 
    shifted = GENOUTS.shifted;
    s_name = [new_dir_name 'shifted_' base_name run_name];
    save(s_name,'shifted'); clear shifted
  end
  if SIMOPTS.save_parents, 
    parents = GENOUTS.parents;
    par_name = [new_dir_name 'parents_' base_name run_name];
    save(par_name,'parents'); clear parents
  end
  if SIMOPTS.save_kills, 
    kills = GENOUTS.kills;
    kills = kills(1:len_p,:);
    kill_name = [new_dir_name 'kills_' base_name run_name];
    save(kill_name,'kills');  clear kills
    if SIMOPTS.save_rivalries, 
      rivalries = GENOUTS.rivalries;
      rvl_name = [new_dir_name 'rivalries_' base_name run_name];
      save(rvl_name,'rivalries'); clear rivalries
    end
  end
else, 
  pop_name = [new_dir_name 'population_' base_name run_name];
  load(pop_name);
  population = [population GENOUTS.population];
  save(pop_name,'population');  clear population
  if SIMOPTS.record,  
    if SIMOPTS.exp_type~=0, 
      tn_name = [new_dir_name 'trace_noise_' base_name run_name];
      load(tn_name);
      trace_noise = [trace_noise; GENOUTS.trace_noise];
      save(tn_name,'trace_noise');  clear trace_noise
    end
    tx_name = [new_dir_name 'trace_x_' base_name run_name];
    load(tx_name);
    trace_x = [trace_x; GENOUTS.trace_x];
    save(tx_name,'trace_x');  clear trace_x
    ty_name = [new_dir_name 'trace_y_' base_name run_name];
    load(ty_name);
    trace_y = [trace_y; GENOUTS.trace_y];
    save(ty_name,'trace_y');  clear trace_y
  end
  if SIMOPTS.save_seeds,  
    tcs_name = [new_dir_name 'trace_cluster_seed_' clus_name run_name];
    load(tcs_name);
    trace_cluster_seed = [trace_cluster_seed; GENOUTS.trace_cluster_seed];
    save(tcs_name,'trace_cluster_seed');  clear trace_cluster_seed
    sd_name = [new_dir_name 'seed_distances_' clus_name run_name];
    seed_distances = [seed_distances; GENOUTS.seed_distances];
    save(sd_name,'seed_distances');
  end
  if size(unique(GENOUTS.land),1)~=1, 
    l_name = [new_dir_name 'land_' base_name run_name];
    load(l_name);
    land = [l GENOUTS.land];
    save(l_name,'land');  clear land
  end
  if size(GENOUTS.shifted,1)~=0, 
    s_name = [new_dir_name 'shifted_' base_name run_name];
    load(s_name);
    shifted = [GENOUTS.shifted shifted];
    save(s_name,'shifted'); clear shifted
  end
  if SIMOPTS.save_parents, 
    par_name = [new_dir_name 'parents_' base_name run_name];
    load(par_name);
    parents = [parents; GENOUTS.parents];
    save(par_name,'parents'); clear parents
  end
  if GENOUTS.save_kills, 
    kill_name = [new_dir_name 'kills_' base_name run_name];
    load(kill_name);
    kills = [kills; GENOUTS.kills];
    save(kill_name,'kills');  clear kills
  end
end
end