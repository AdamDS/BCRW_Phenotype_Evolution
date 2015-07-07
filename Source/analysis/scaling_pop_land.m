%% scaling_pop_land.m *******************************************************

global SIMOPTS;
m = 0; c = 'ybrkybrk'; c = [c c c c c];
op = critical_points(1);
dm = critical_points(2);
mu = critical_points(3);
SIMOPTS.op = op;
SIMOPTS.dm = dm;
SIMOPTS.mu = mu;
lbmslS = zeros(length(basic_map_sizes),length(SIMS));
avg_lps = lbmslS;
avg_end_lps = lbmslS;
avg_dps = lbmslS;
avg_end_dps = lbmslS;

lbms = lbmslS(:,1);
AVG_LPS = lbms;
STD_LPS = lbms;
MIN_LPS = lbms;
MAX_LPS = lbms;
RAT_LPS = lbms;
SRAT_LPS = lbms;
AVG_NGENS_LPS = lbms;
STD_NGENS_LPS = lbms;
AVG_END_LPS = lbms;
STD_END_LPS = lbms;

lbmsN = zeros(length(basic_map_sizes),NGEN);
AVG_GEN_LPS = lbmsN;
STD_GEN_LPS = lbmsN;

num_sims = zeros(length(basic_map_sizes),1);
tic;
for bms = basic_map_sizes,  
  m = m +1;
  SIMOPTS.basic_map_size = [bms bms]; 
  if default_density, SIMOPTS.IPOP = make_IPOP(SIMOPTS.basic_map_size); end
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir,do_cd);
  gen_pops = zeros(length(SIMS),NGEN);
  for run = SIMS, 
    run_name = int2str(run);
    pop_name = ['population_' base_name run_name];
    [p,go] = try_catch_load(pop_name,1,1);
    if go,  
      population = p.population;  clear p
      this_NGEN = length(find(population));
      g = find(population>=limit);
      avg_lps(m,run) = mean(population);
      if do_end_plot, if this_NGEN>transience,  
        avg_end_lps(m,run) = mean(population(transience+1:end));  
      else, 
        avg_end_lps(m,run) = 0;
      end,  end
    end
  end %for run
  exists = find(avg_lps(m,:));
  num_sims(m) = length(exists);
  if exists,  
    AVG_LPS(m) = mean(avg_lps(m,exists));
    STD_LPS(m) = std(avg_lps(m,exists));
    AVG_END_LPS(m) = mean(avg_end_lps(m,exists));
    STD_END_LPS(m) = std(avg_end_lps(m,exists));
  end
  
  [area(m),xl(m),yl(m)] = landscape_measures(SIMOPTS.basic_map_size);
  AVG_LDS(m) = AVG_LPS(m)./area(m);
  STD_LDS(m) = STD_LPS(m)./area(m);
  AVG_END_LDS(m) = AVG_END_LPS(m)./area(m);
  STD_END_LDS(m) = STD_END_LPS(m)./area(m);
end %for bms
toc;

if do_pop_land_plot,  
  figure(185600000+length(basic_map_sizes)+(1000*mu));
  subplot(2,2,1);
  plot(area,AVG_LPS,'-x');
  subplot(2,2,3);
  plot(area,STD_LPS,'-x');
  subplot(2,2,2);
  plot(area,AVG_END_LPS,'-x');
  subplot(2,2,4);
  plot(area,STD_END_LPS,'-x');
end
if do_den_land_plot,  
  figure(185700000+length(basic_map_sizes)+(1000*mu));
  subplot(2,2,1);
  plot(area,AVG_LDS,'-x');
  subplot(2,2,3);
  plot(area,STD_LDS,'-x');
  subplot(2,2,2);
  plot(area,AVG_END_LDS,'-x');
  subplot(2,2,4);
  plot(area,STD_END_LDS,'-x');
end

lx = log10(xl); lALPS = log10(AVG_LPS);
lAELPS = log10(AVG_END_LPS);
[eslope,einter,esigm,esigb] = linear_fit(lx,lAELPS',[],1);
% [slope,inter,sigm,sigb] = linear_fit(lx,lALPS',[],1);