global SIMOPTS;
i = 0;
AVG_LCRAT = zeros(length(mutability),1);
STD_LCRAT = zeros(length(mutability),1);
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  mlcrat = zeros(length(SIMS),1);
  i = i +1;
  for run = SIMS
    run_name = int2str(run);
    if limit~=3
      load(['population_' base_name run_name]);
      load(['num_clusters_' base_name int2str(limit) '_limit_' run_name]);
      load(['orgsnclusters_' base_name int2str(limit) '_limit_' run_name]);
    else
      load(['population_' base_name run_name]);
      load(['num_clusters_' base_name run_name]);
      load(['orgsnclusters_' base_name run_name]);
    end
    u = 0;  v = 0;  nu = 0; nv = 0;
    lcrat = zeros(length(population),1);
    for gen = 1:length(population)
      u = v +1;
      v = sum(population(1:gen));
      nu = nv +1;
      nv = sum(num_clusters(1:gen));
      pop = population(gen);
      nc = num_clusters(gen);
      onc = orgsnclusters(nu:nv);
      monc = max(onc);
      lcrat(gen) = monc/pop; %ratio of largest cluster to population
    end
    if mu==7.0 && run==1
      figure, plot(lcrat);
    end
    mlcrat(run) = mean(lcrat(find(lcrat)));  lcrat = [];
  end
  AVG_LCRAT(i) = mean(mlcrat);
  STD_LCRAT(i) = std(mlcrat); mlcrat = [];
end
end
end
if do_lcrat_plot==1
  figure(5300 + limit); errorbar(mutability,AVG_LCRAT,STD_LCRAT);
  title(make_title_name(generalize_base_name(base_name),''));
  xlabel('\mu');  ylabel('ratio of largest cluster size to population');
end