%% main_analysis.m *******************************************************
% This function is based on the main_babies template which means setting up
% analysis for a particular simulation may be done by choosing the typical
% options. Note that at the end, there are options for gathering population
% or num_clusters data. More updates which will allow for other analysis
% options may be added there as well.


% pctconfig('portrange',[27370 27470])
%% START CLEAN & RANDOM
clear; clc; close all;
rng('shuffle','twister');
% POLICY = system_dependent('RemoteCWDPolicy', 'None');

%% PARAMETERS
% Data creation options
source = 'G:\Babies_Root\';%'G:\Babies_Root\';%Bortas\Bacterial\Mu_x\Uniform\2_Flatscape\Data';
  do_cd = 0;
output = 'C:\Users\amviot\Desktop\';
for_external = 0;
OS = 0;
write_over = 1;
pool_size = 0;
no_bio = 0;
make_dir = 0;
split = 10;

% What data to create
only_lt = 0; %record only times to fixation if 1 and disregard all other data
only_relax = 0;

% Simulation numbers
SIMS = 1:10;%1:100;%simulation identifiers used in simulation looping
SIMrange = [int2str(SIMS(1)) '_' int2str(SIMS(end))];

% Initial parameter settings
NGEN = 2000;%10^5;%2*10^3; %number of generations to run
transience = 500;
IPOP = 300;%30000; %number of initial population
  do_default_density = 1;
limit = 3; %minimum cluster size/extinction population +1
loaded = 0; %choose to load predefined variables, 0 = no load, 1 = load load_name variables
load_name = ['']; %string which identifies predefined variables babies and basic_map
append = 0;

% Birth settings
  % Reproduction options
  distribution = 0; %offspring distrubtion: 0=uniform, 1=normal
  reproduction = 1; %assortative mating = 0, bacterial cleaving = 1, random mating = 2
  % Mutability option
  exp_type = 0; %same mutability = 0; competition = 1; duel = 2
  % Single Mutability settings
  dmu = 0.01;  mutability = [0.30:dmu:0.36];%0.30:dmu:0.36;%[0.25:0.01:0.33 0.331:dmu:0.34 0.36:0.01:0.45];
  relaxed = [0];% 0 1200 1800 2500 2700 2900];
  % Two mu competition settings
  bi = [1.31 2.31; 150 150]; %[mu1 mu2; IPOP1 IPOP2];
  % Mutability competition settings
  range = 7; % 0 to range possible mutabilities

% Death settings
  % Local options
  dop = 0.25; overpop = 0.25; %if closer than this distance, overpopulated, and baby dies
  random_walk = 0; %0 = coalescing, 1 = annihilating
  % Global options
  ddm = 0.1;  death_max = 0.7; %percent of random babies dying varies from 0 to this value.
  indiv_death = 0;  %random percentage of entire pop dies = 0; individual probability = 1

% Landscape options
  % Heights
  shock = 0; %Set this to 1 to generate a new random map every landscape_movement generations
    shock_heights = [1 1]; %min and max landscape heights during shock
    shock_duration = 2; %duration of a shock (generations)
    shock_repeat = 0; %repeated shocks, 0 = one shock, 1 = shocks every landscape_movement generations
    shock_over = 0; %flag for when single shock has occurred, 0 = hasn't finished, 1 = finished shock
  landscape_heights = [2 2]; %min and max of landscape; for flat landscapes, only min is taken
  %If both values are the same, then a flatscape will be generated.
  % Movement
  landscape_movement = 2*NGEN; %land moves every "landscape_movement" generations
    %There is only a default shift of 1 basic_map row, we could add this in if 
    %we really want to, but it may be unnecessary.
  % Size
  basic_map_sizes = 12;%[8 10 12]; %X and Y lengths for the basic map size
  basic_map_size = [12 12];
  linear = 0;
    %May just put in a check on whether one basic_map_size values is 1
  periodic = [0 0]; %periodic boundary conditions for the x & y coords, respectively
    %if [1 0] or [0 1], then cylindrical boundaries wrapping x or y edges, respectively
    %if [1 1], then the landscape becomes toroidal
  %does not yet work with varied landscape iterpolations
  INFRAT = 0.989;
  
  if do_default_density,  defIPOP = IPOP; end
  
  global SIMOPTS;
  SIMOPTS = struct('source',source,'OS',OS,'pool_size',pool_size,'only_lt',only_lt,...
    'no_bio',no_bio,'write_over',write_over,'split',split,'for_external',for_external,...
    'SIMS',SIMS,'NGEN',NGEN,'IPOP',IPOP,'defIPOP',defIPOP,'limit',limit,...
    'loaded',loaded,'load_name',load_name,...
    'op',[],'random_walk',random_walk,'dm',[],'indiv_death',indiv_death,...
    'shock',shock,'shock_heights',shock_heights,'shock_duration',shock_duration,...
    'shock_repeat',shock_repeat,'shock_over',shock_over,...
    'landscape_movement',landscape_movement,'landscape_heights',landscape_heights,...
    'basic_map_sizes',basic_map_sizes,'basic_map_size',basic_map_size,'linear',linear,...
    'periodic',periodic,'distribution',distribution,'reproduction',reproduction,...
    'exp_type',exp_type,'mu',[],'bi',bi,'range',range);
    
if pool_size>0 && matlabpool('size')==0, matlabpool open; end
total_params = length(mutability)*length(overpop)*length(death_max);

%% What to check
             do_check_data_set = 0;   do_passes = 0;  do_fails = 1;   do_DNE = 0;   do_only_DNE = 0;
                                      do_pause_each = 0;  do_pause_set = 0;
          do_check_duplication = 0;
       do_check_seed_distances = 0;                   try_sqrt = 0;
          
%% What to fix
      do_fix_data = 0;    do_fix_population = 1;
%% What to record
record_data = 0;
export_coords = 0;
record_abundances = 0;
record_R = 0;
record_rank_abundance = 0;
record_kills_v_pop = 0;

%% What to analyze
        %Data & recording                 %Visualizations                    %Options
                do_populations = 0;        do_populations_plot = 0;             by_generation = 1;
                                               do_density_plot = 0;                 by_window = 0;
                                                   do_end_plot = 0;                     END = 100;
                                                   do_std_plot = 0;
                                             do_landscape_plot = 0;
               do_num_clusters = 0;       do_num_clusters_plot = 1;
           do_scaling_pop_land = 0;           do_pop_land_plot = 1;       
                                              do_den_land_plot = 0;           critical_points = [0.25,...
                                                                                          0.70,0.337];
              do_connectedness = 1;     do_connectedness_plots = 1;
             do_seed_distances = 0;  
           do_Clark_Evans_test = 0;        do_population_level = 1;
                                              do_cluster_level = 0;
                                     do_nearest_neighbors_plot = 1;              
                                        do_parent_distances_3d = 0; 
                                                  do_dist_plot = 0;              
                                                     do_R_plot = 1;                      ddpr = 1;
                                                                             do_modified_area = 0;
                 do_abundances = 0;         do_abundances_plot = 0;      do_abundances_single = 0;%  
                                             do_abundances_log = 0;     do_abundances_general = 1;
                                             do_abundances_avg = 0;       
                                              do_abundances_3d = 0; 
                                       do_abundances_log_error = 0;   
            do_rank_abundances = 0;
                      do_kills = 0;              do_kills_plot = 1;                 all_kills = 1;
                do_kills_v_pop = 0;       do_kills_v_pop_plots = 1;
              do_rivalry_v_pop = 0;
                do_plot_indivs = 0;                 make_video = 0;          plot_mu = mutability;%
                                                                               plot_run = SIMS(1);
                                                                plot_gens = [1 ceil(NGEN/2) NGEN];
                                                                                     subplots = 1;
                                                                                subplotrc = [3 3];
                                                                          generate_paper_plot = 0;
                                                                            backgroundcolor = 'w';
                                                                        plot_largest_clusters = 1;
                                                                             color_code = [1 2 3];
                 do_hist_coord = 0;                                       plot_hc_mu = mutability;%
                                                                                  plot_hc_run = 1;
                                                                                plot_hc_gen = 100;
                do_corr_pop_nc = 0;               do_corr_plot = 1;
                      do_lcrat = 0;              do_lcrat_plot = 0;
                         do_cp = 0;                 do_cp_plot = 0;
              do_indiv_lineage = 0;        do_descendants_mesh = 0;
      do_indiv_cluster_lineage = 0;
            do_cluster_lineage = 0;                do_cl_hists = 0;
                                                   do_cl_stats = 0;          do_cl_stat_plots = 0;
                                                            to = 1;
                                                          from = 1;
   do_plot_lineage_proportions = 0;
                do_percolation = 0;        do_percolation_plot = 1;      percolation_limit = 0.95;
                       do_beta = 0;
% do_cluster_number_distribution = 0;
%              do_gyration_radii = 0;
%            do_mass_frequencies = 0;
%             do_cluster_numbers = 0;
%           do_correlation_lengths = 0;               
do_correlation_lengths_distribution = 0;  do_correlation_lengths_plot = 1;
       do_infinite_probability = 0;
         do_cluster_spans_plot = 0;
%                 do_nu_parallel = 0; mucrit = 0.335; delta = 0.462; %0.462 +- 0.033 from mu=0.335 for 100sims
%                                     nu_parallel = 1.49; %0.335 1.05
%                                     window = 10;
%                                     t = 1:NGEN;
%                                     use = 1:10:NGEN;
         do_critical_exponents = 0;           do_plot_collapse = 0;              offset = 0;   
                                                    do_plot_Pb = 1;                  dt = 1;  
                                                                                 eta = 0.05;
                                                                              amM = [10 260];
                                                                             fTmM = [10 510];   
                                                                             fSmM = [10 510];
% do_spatial_correlation_lengths = 0;
%     do_correlation_length_rate = 0;
%      do_characteristic_lengths = 0;do_characteristic_diameters = 1;
%                                         do_correlation_lengths = 1;                N_samples = 10;
%                                                                                 epsilon = overpop;
%                                                                                    de = overpop*2;
%                                                                                      EPSILON = 45;
  do_cluster_density_dimension = 0;
          do_cluster_dimension = 0;                 samples = 1000;
                                                begin_of_end = 200;
                                       cluster_dimension_IR = 0.98;
            do_fisher_exponent = 0;                                            do_fisher_plot = 1;
                do_filled_area = 0;
                    do_degrees = 0;            do_degrees_plot = 0;
                do_path_length = 0;           num_final_gens = 100;
      do_lifetimes_to_survival = 0;                              do_survival_probability_plot = 1;
   do_transition_probabilities = 0;                   binsize = 25;        do_transition_mesh = 0;
                                                                     do_transition_ratio_plot = 1;
               do_plot_lineage = 0;

%% The Analysis
tic
if only_lt, 
  write_lifetimes;
else
  if do_populations, get_populations;  end %updated 11/1/12

  if do_num_clusters, get_num_clusters;  end %updated 11/1/12

  if record_data,
    if do_populations, write_populations;  end
    if do_num_clusters, write_clusters; end
    if do_populations && do_num_clusters, write_lifetimes; end
  end   %updated 11/1/12

  if do_scaling_pop_land, scaling_pop_land; end
  
  if do_seed_distances, make_seed_distances;  end
  
  if do_Clark_Evans_test, 
    get_R;
    if record_R, write_R;  end
  end  %updated 11/1/12
  
  if do_abundances, get_abundances_curve;  end

  if do_rank_abundances, get_rank_abundance;  end

  if do_corr_pop_nc, corr_pop_nc;  end
  
  if do_plot_indivs, plot_indivs;  end
  
  if do_hist_coord, hist_coord;  end
  
  if do_kills, get_kills;  end

  if do_kills_v_pop, kills_v_pop;  end
  
  if do_rivalry_v_pop, rivalry_v_pop;  end

  if do_degrees, get_degrees;  end
  
  if do_connectedness, get_connectedness;  end   %updated 2/22/12
  
  if do_indiv_lineage, get_indiv_lineage;  end
  
  if do_indiv_cluster_lineage, get_indiv_cluster_lineage;  end

  if do_cluster_lineage, 
    if to, build_Cluster_Lineage_to; end
    if from, build_Cluster_Lineage_from; end
    if do_cl_hists,  get_cluster_lineage_hists;  end
    if do_cl_stats, get_lineage_stats;  end
  end
  
  if do_plot_lineage, plot_lineage; end
  
  if do_plot_lineage_proportions, plot_lineage_proportions; end
  
  if do_correlation_lengths_distribution, get_correlation_lengths_distribution; end
  
  if do_infinite_probability, get_infinite_probability; end
  
%   if do_nu_parallel,  get_nu_parallel;  end

  if do_critical_exponents, get_critical_exponents; end
  
  if do_cluster_spans_plot, get_cluster_spans_plot; end
  
  if do_fisher_exponent,  get_fisher_exponent;  end
   
  if do_path_length, get_path_length; end
  
  if do_lifetimes_to_survival,  get_lifetimes_to_survival_probabilities; end
  
  if do_transition_probabilities, get_transition_probabilities; end
  
  if do_percolation, percolation_lengths;  end
  
  if do_beta, get_beta; end
%%
  if do_check_seed_distances,  check_seed_distances; end
  if do_check_data_set,  check_data_set; end
  if do_check_duplication,  check_duplication;  end
%%
  if do_fix_data, fix_data; end
end
toc
if pool_size>0, matlabpool close; end

%% Old code
%   
%   if do_gyration_radii, get_gyration_radii; end
% 
%   if do_cluster_number_distribution,  make_cluster_number_distribution; end
%   
%   if do_gyration_radii, get_gyration_radii; end
% 
%   if do_mass_frequencies, build_mass_frequencies; end
%   
%   if do_cluster_numbers,  build_cluster_numbers;  end
% 
%   if do_correlation_lengths,  build_correlation_lengths;  end
%
%   if do_spatial_correlation_lengths, get_spatial_correlation_lengths;  end
%   
%   if do_characteristic_lengths,  get_characteristic_length_distributions;  end
%   
%   if do_lcrat, get_ratio_of_largest_cluster; end
% 
%   if do_cp, continuum_percolation; end
%   
%   if do_cluster_density_dimension, get_cluster_density_dimension;  end
%     
%   if do_filled_area, get_filled_area;  end
%   
%   if do_correlation_length_rate, get_correlation_length_growth_rate; end