%% get_beta.m *******************************************************
% critical_points = [0.25 0.7 0.334];
global SIMOPTS;
m = 0; c = 'ybrkybrk'; c = [c c c c c];

lmlS = zeros(length(mutability),length(SIMS));
avg_pops = lmlS;
avg_end = lmlS;

lm = lmlS(:,1);
AVG_POPS = lm;
STD_POPS = lm;
AVG_END = lm;
STD_END = lm;

MU = critical_points(3) -mutability';
num_sims = zeros(length(mutability),1);
tic;
% for bms = basic_map_sizes,  
SIMOPTS.basic_map_size = basic_map_size; 
if default_density, SIMOPTS.IPOP = make_IPOP(SIMOPTS.basic_map_size); end
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability,  SIMOPTS.mu = mu;
  m = m +1;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir,do_cd);
  r = 0;
%   gen_pops = zeros(length(SIMS),NGEN);
  for run = SIMS, 
    r = r +1;
    run_name = int2str(run);
    pop_name = ['population_' base_name run_name];
    [p,go] = try_catch_load(pop_name,1,1);
    if go,  
      population = p.population;  clear p
      this_NGEN = length(find(population));
      g = find(population>=limit);
      avg_pops(m,r) = mean(population);
      if do_end_plot, if this_NGEN>transience,  
        avg_end(m,r) = mean(population(transience+1:end));  
      else, 
        avg_end(m,r) = 0;
      end,  end
    end
  end %for run
  exists = find(avg_pops(m,:));
  num_sims(m) = length(exists);
  if exists,  
    AVG_POPS(m) = mean(avg_pops(m,exists));
    STD_POPS(m) = std(avg_pops(m,exists));
    AVG_END(m) = mean(avg_end(m,exists));
    STD_END(m) = std(avg_end(m,exists));
  end
end %for mu
end %for dm
end %for op
[area,xl,yl] = landscape_measures(SIMOPTS.basic_map_size);
AVG_DENS = AVG_POPS./area;
STD_DENS = STD_POPS./area;
AVG_END_DENS = AVG_END./area;
STD_END_DENS = STD_END./area;
toc;

lmu = log10(MU);
lden = log10(AVG_DENS);
[beta,inter,sigm,sigb] = linear_fit(lmu,lden,[],1);