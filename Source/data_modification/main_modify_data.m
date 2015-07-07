%% modify_data.m *******************************************************
% This may be used to modify various data. The template is that of main_babies.m and
% main_clusters.m, so you can look to those help comments for setting up your
% modifications.

%% START CLEAN
clear; clc; %close all;

%% PARAMETERS

% Data creation options
  opt = 0; %change directories = 0, change name only = 1
  source = 'G:\Babies_Root\';%'C:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\10_NGEN\Data\';
  do_cd = 0;
  write_over = 0;
  OS = 0;
  make_dir = 0;
  do_cd = 0;
  split = 0;

  % What to do
  no_bio = 0;
  only_lt = 0; %record only times to fixation if 1 and disregard all other data
  only_relax = 0;
  do_simulations = 1; %generates population, trace_x, trace_y, trace_noise, trace_cluster_seed*
  do_clustering = 1;  %generates trace_cluster_seed* & allows building and locating clusters
    do_build_clusters = 1;  %generates num_clusters, trace_cluster, orgsnclusters, 
    do_locate_clusters = 1; %generates centroid_x, centroid_y, cluster_diversity
    do_build_diameters = 1;
  do_genealogies = 0; %generates lineage info of original population & of the clusters
    do_indiv_lineage = 1; %generates num_descendants
    do_indiv_cluster_lineage = 1; %generates num_descendant_clusters, descendant_clusters
    do_cluster_lineage = 1; %generates num_clusters_produced, clusters produced, 
                            %num_clusters_fused, clusters_fused

% Simulation numbers
SIMS = [1:1000]; %simulation identifiers used in simulation looping

% Initial parameter settings
NGEN = 2*10^3; %number of generations to run
IPOP = 300; %number of initial population
  default_density = 0;
relax = 0;
limit = 3; %minimum cluster size/extinction population +1
loaded = 0; %choose to load predefined variables, 0 = no load, 1 = load load_name variables
load_name = ['']; %string which identifies predefined variables babies and basic_map

% Death settings
  % Local options
  overpop = 0.25; %if closer than this distance, overpopulated, and baby dies
  random_walk = 0; %0 = coalescing, 1 = annihilating
  % Global options
  ddm = 0;  death_max = [0.7]; %percent of random babies dying varies from 0 to this value.
  indiv_death = 0;  %random percentage of entire pop dies = 0; individual probability = 1

% Landscape options
shock = 0; %Set this to 1 to generate a new random map every landscape_movement generations
shock_heights = [1 1]; %min and max landscape heights during shock
shock_duration = 5; %duration of a shock (generations)
shock_repeat = 0; %repeated shocks, 0 = one shock, 1 = shocks every landscape_movement generations
shock_over = 0; %flag for when single shock has occurred, 0 = hasn't finished, 1 = finished shock
landscape_movement = 200; %land moves every "landscape_movement" generations
  %There is only a default shift of 1 basic_map row, we could add this in if 
  %we really want to, but it may be unnecessary.
landscape_heights = [2 2]; %min and max of landscape; for flat landscapes, only min is taken
  %If both values are the same, then a flatscape will be generated.
basic_map_size = [12 12]; %X and Y lengths for the basic map size
linear = 0;
  %May just put in a check on whether one basic_map_size values is 1
periodic = [0 0];

% Reproduction options
distribution = 0; %offspring distrubtion: 0=uniform, 1=normal
reproduction = 1; %assortative mating = 0, bacterial cleaving = 1, random mating = 2

% Mutability option
exp_type = 0; %same mutability = 0; competition = 1; duel = 2

% Single Mutability settings
dmu = 0.001;  mutability =[0.250:dmu:0.320];
relaxed = [1200];% 0 1200 1800 2500 2700 2900];

% Two mu competition settings
bi = [1.31 2.31; 150 150]; %[mu1 mu2; IPOP1 IPOP2];

% Mutability competition settings
range = 7; % 0 to range possible mutabilities

global SIMOPTS;
tic;
for op = overpop
for dm = death_max
for mu = mutability
  bad_bn = int2str(mu*100);
  %% This has to be exact and set by you! The underscore at the end is necessary, so don't forget it!
%   bad_dir = ['C:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\10_NGEN\Data\'];
  
  % Bundle simulation options into SIMOPTS
  SIMOPTS = struct('source',source,'OS',OS,'split',split,...
    'only_lt',only_lt,'only_relax',only_relax,'relaxed',relaxed,...
    'do_simulations',do_simulations,'do_clustering',do_clustering,...
    'do_build_clusters',do_build_clusters,'do_locate_clusters',do_locate_clusters,...
    'do_build_diameters',do_build_diameters,...
    'do_genealogies',do_genealogies,'do_indiv_lineage',do_indiv_lineage,...
    'do_cluster_lineage',do_cluster_lineage,'do_indiv_cluster_lineage',do_indiv_cluster_lineage,...
    'SIMS',SIMS,'NGEN',NGEN,'IPOP',IPOP,'limit',limit,'loaded',loaded,'load_name',load_name,...
    'op',op,'random_walk',random_walk,'dm',dm,'indiv_death',indiv_death,...
    'shock',shock,'shock_heights',shock_heights,'shock_duration',shock_duration,...
    'shock_repeat',shock_repeat,'shock_over',shock_over,...
    'landscape_movement',landscape_movement,'landscape_heights',landscape_heights,...
    'basic_map_size',basic_map_size,'linear',linear,'periodic',periodic,...
    'distribution',distribution,'reproduction',reproduction,'exp_type',exp_type,'mu',mu,...
    'bi',bi,'range',range);
  %build the simulation names for which you have chosen to generate cluster data
  [base_name,bad_dir] = NameAndCD(make_dir,do_cd);
  for run = SIMS, 
    run_name = int2str(run);
    bad_name = [bad_dir base_name]; %splits this if opt==0
    fprintf(['changing directory of %s \n'],bad_dir);
    good_dir = split_cd(bad_dir,run,10,1,0);
    fprintf(['to %s \n'],good_dir);
    good_name = [good_dir base_name]; %splits this if opt==0
%     fprintf(['changing name of ' bad_name run_name '\n']);
    %% Include any data modification functions here
    movefile_data_set(bad_name,good_name,run_name,opt); %
    %     fix_population(base_name,run_name,limit);
%     new_name = 'a'; %test/debug name
%     if exist(['nn_dist_' bad_name run_name '.mat'])==2
%       load(['nn_dist_' bad_name run_name]);  
%       nn_dist = sqrt(nn_dist); 
%       save(['nn_dist_' base_name run_name '.mat'],'nn_dist');
%     end
%     if exist(['mnd_' bad_name run_name '.mat'])==2
%       load(['mnd_' bad_name run_name]);  
%       mnd = sqrt(mnd); 
%       save(['mnd_' base_name run_name '.mat'],'mnd');
%     end
%     if exist(['trace_cluster_seed_' bad_name run_name '.mat'])==2
%       movefile(['trace_cluster_seed_' bad_name run_name '.mat'],['trace_cluster_seed_' base_name run_name '.mat']);
%     end
%     delete_data_set(bad_name,run_name);
%     rename_data_set(bad_name,base_name,run_name);
% Uncommenting delete_data_set will get rid of the "badly named" file set, leaving your fresh, new data set alone.
%     delete_data_set(bad_name,run_name);
  end
end
end
end
toc;