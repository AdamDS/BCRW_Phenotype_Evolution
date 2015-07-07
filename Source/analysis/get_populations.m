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
m = 0; c = 'ybrkybrk'; c = [c c c c c];
for bms = basic_map_sizes,  
  m = m +1;
  SIMOPTS.basic_map_size = [bms bms]; 
  
i = 1;

lmlS = zeros(length(mutability),length(SIMS));
avg_pop = lmlS;
ngens = lmlS;

lm = lmlS(:,1);
AVG_POPS = lm;
STD_POPS = lm;
AVG_NGENS = lm;
STD_NGENS = lm;

lmN = zeros(length(mutability),NGEN);
AVG_GEN_POPS = lmN;
STD_GEN_POPS = lmN;

num_sims = zeros(length(mutability),1);
sim_exists = lmlS;

clear lmlS lm lmN

i = 0;
tic;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  [base_name,dir_name] = NameAndCD(make_dir,do_cd);
  i = i +1;
  r = 0;
  if only_lt, 
    [p,go] = try_catch_load([dir_name 'populations_' base_name SIMrange],1,1);
    if go,  [lt,go] = try_catch_load([dir_name 'lifetimes_' base_name SIMrange],1,1);
    if go,  
      populations = p.populations;  clear p
      lifetimes = lt.lifetimes; clear lt
      gens = [ones(lenght(lifetimes),1) lifetimes];
      AVG_POPS(i) = mean(;
      STD_POPS(i) = std(populations);
      AVG_GEN_POPS(i,:) = mean(populations);
      STD_GEN_POPS(i,:) = std(populations);
      ngens(i,:) = lifetimes;
    end
    end
  else, 
    %build the simulation names for which you have chosen
    gen_pops = zeros(length(SIMS),NGEN);
    for run = SIMS, 
      r = r +1;
      new_dir_name = split_cd(dir_name,run,split,make_dir,do_cd);
      run_name = int2str(run);
      exp_name = [base_name run_name];
      pop_name = [new_dir_name 'population_' exp_name];
      [p,go] = try_catch_load(pop_name,1,1);
      if go, 
        num_sims(i) = num_sims(i) +1;
        sim_exists(i,r) = 1; %track which simulation data is used
        population = p.population;  clear p
        if ~loaded || no_bio, 
          g = find(population>=limit);
          actual_ngen = max(g);
        else, 
          g = find(population);
          actual_ngen = length(g);
        end
        ngens(i,r) = actual_ngen;
        avg_pop(i,run) = mean(population(g));
        gen_pops(run,1:length(g)) = population(g);
      end
    end
  end
  use = find(sim_exists(i,:)); these = 1:length(use);
  gen_pops(~use,:) = [];

  [rgood(i),cgood] = size(gen_pops);

  AVG_POPS(i) = mean(avg_pop(i,:));
  STD_POPS(i) = std(avg_pop(i,:));

  AVG_GEN_POPS(i,:) = mean(gen_pops,1);
  STD_GEN_POPS(i,:) = std(gen_pops);
end
end
end
toc;
[A(m),xl(m),yl(m)] = landscape_measures(SIMOPTS.basic_map_size);
area = A(m);
fn = 0; lbms(m) = bms;

%% Plot populations by control parameter
if do_populations_plot, 
  fn1 = 1;  
  plot_order_v_control(mutability,AVG_POPS,STD_POPS,fn1,...
    '\mu','populations',c(m),base_name);  
  figure(fn1);  hold on;  if do_std_plot, figure(fn1+1);  hold on;  end
end
%% Plot populations by generation
if do_populations_plot==1 && by_generation==1, 
  figure(3);
  tn = make_title_name(generalize_base_name(base_name),'');
  plot([1:NGEN],AVG_GEN_POPS);  title(tn,'FontSize',16);  hold on;
  xlabel('generation','FontSize',14); ylabel('<population>','FontSize',14);
  xlim([1 NGEN]); 
end
end
if do_populations_plot, 
  figure(fn1);  legend(int2str(basic_map_sizes'));
  if do_std_plot, figure(fn1+1);  legend(int2str(basic_map_sizes'));  end
end
%% Plot landscape by mutability
dis = 3;
if do_landscape_plot,  
  if do_populations_plot, 
    fn5 = 10;  %if do_std_plot, fn5(2) = 100;  else, fn5(2) = []; end
    plot_order_v_control(xl,land_pops(dis,:),land_std_pops(3,:),2000+fn5,...
      'land','populations',c(m),base_name);
  end
end