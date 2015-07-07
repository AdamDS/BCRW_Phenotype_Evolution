global SIMOPTS;
i = 0;
AVG_DIST_RO_L = zeros(length(mutability),1);
STD_DIST_RO_L = zeros(length(mutability),1);
AVG_POP_DENS = zeros(length(mutability),1);
STD_POP_DENS = zeros(length(mutability),1);
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  mpopd = zeros(length(SIMS),1);
  i = i +1;
  for run = SIMS
    run_name = int2str(run);
    if limit~=3
      load(['population_' base_name run_name]);
%       load(['trace_cluster_seed_' base_name int2str(limit) '_limit_' run_name]);
      load(['trace_cluster_seed_' base_name run_name]);
      load(['trace_x_' base_name run_name]);
      load(['trace_y_' base_name run_name]);
    else
      load(['population_' base_name run_name]);
      load(['trace_cluster_seed_' base_name run_name]);
      load(['trace_x_' base_name run_name]);
      load(['trace_y_' base_name run_name]);
    end
    u = 0;  v = 0;
    mdrol = zeros(length(population),1);
    bmsx = basic_map_size(1);  bmsy = basic_map_size(2);
    lsx = (((bmsx*2)-1)*2)-1; lsy = (((bmsy*2)-1)*2)-1;
    for gen = 1:length(population)
      u = v +1; 
      v = sum(population(1:gen));
      ltcs = trace_cluster_seed(u:v,limit-2);
      x = trace_x(u:v); y = trace_y(u:v);
      dist_RO_L = ((x -x(ltcs)).^2 +(y -y(ltcs)).^2).^0.5; %distance between RO and least alike seed member
      mdrol(gen) = mean(dist_RO_L);
    end
    pop_dens = population /(lsx *lsy);
    mmdrol(run) = mean(mdrol);
    mpopd(run) = mean(pop_dens);
  end
  AVG_DIST_RO_L(i) = mean(mmdrol);
  STD_DIST_RO_L(i) = std(mmdrol);
  AVG_POP_DENS(i) = mean(mpopd);
  STD_POP_DENS(i) = std(mpopd);
end
end
end
if do_cp_plot==1
  figure(5120000 + limit); 
  subplot(3,1,1);
  plot(AVG_POP_DENS,AVG_DIST_RO_L,'x');
  xlabel('<\rho_{pop}>'); ylabel('<\Delta_{RO,L}>');
  title(make_title_name(generalize_base_name(base_name),''));
%   figure(5130000 +limit);
  subplot(3,1,2);
  errorbar(mutability,AVG_POP_DENS,STD_POP_DENS);
  xlabel('\mu');  ylabel('<\rho_{pop}>');
%   figure(5140000 +limit);
  subplot(3,1,3);
  errorbar(mutability,AVG_DIST_RO_L,STD_DIST_RO_L);
  xlabel('\mu');  ylabel('<\Delta_{RO,L}>');
end