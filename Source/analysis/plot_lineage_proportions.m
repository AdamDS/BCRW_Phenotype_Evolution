%% plot_lineage_proportions.m
%
this_script = 'plot_lineage_proportions';
fprintf([this_script '\n']);
global SIMOPTS;
limit = SIMOPTS.limit;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  [base_name,dir_name] = NameAndCD(0,do_cd);
  
  for run = SIMS
    run_name = int2str(run);
    [p,go] = try_catch_load(['population_' base_name run_name],1,1);
    if go,  [nd,go] = try_catch_load(['num_descendants_' base_name run_name],go,1);
    if go,  
    population = p.population;  clear p
    num_descendants = nd.num_descendants; clear nd
    ipop = SIMOPTS.IPOP;
    ngen = length(find(population>=limit));
    rat = zeros(ipop,ngen);
    rat(:,1) = ones(ipop,1)./population(1);
    for gen = 2:ngen
      rat(:,gen) = num_descendants(:,gen-1)./population(gen);
    end
    spacings = cumsum(rat);
    end
    end
    figure(470000000 +100000*mu +1000*dm +run);
    plot(spacings');  hold on;  plot(1:ngen,zeros(ngen,1)); hold off;
    xlim([1 ngen]); ylim([0 1]);
    xlabel('generation'); ylabel('propotion of lineages');
    title(make_title_name(base_name,run_name));
  end
end
end
end