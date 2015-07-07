%% make_cluster_number_distribution.m
%
%
global SIMOPTS;
this_script = 'make_cluster_number_distribution';
limit = SIMOPTS.limit;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  [base_name,dir_name] = NameAndCD(0,0);
  [clus_name] = cluster_name(base_name);
  NS = [];
  fprintf('Attempting %s for %s \n',this_script,clus_name(1:end-1));
  cnd_name = make_cluster_number_distribution_name(dir_name,clus_name,SIMS);
  if ~mat_exist(cnd_name) || SIMOPTS.write_over,  
  for run = SIMS, 
    run_name = int2str(run);
    
    new_dir_name = split_cd(dir_name,run,split,0,0);
    
    fprintf(' %s of %s \n',[clus_name run_name],int2str(length(SIMS)));
    
    [p,go] = try_catch_load([new_dir_name 'population_' base_name run_name],1,1); 
    if go,  [nc,go] = try_catch_load([new_dir_name 'num_clusters_' clus_name run_name],1,1); 
    if go,  [onc,go] = try_catch_load([new_dir_name 'orgsnclusters_' clus_name run_name],1,1);
    if go,  [tc,go] = try_catch_load([new_dir_name 'trace_cluster_' clus_name run_name],1,1);
    if go,  [tx,go] = try_catch_load([new_dir_name 'trace_x_' base_name run_name],1,1);
    if go,  [ty,go] = try_catch_load([new_dir_name 'trace_y_' base_name run_name],1,1);
    if go,
    population = p.population;  clear p,  
    num_clusters = nc.num_clusters; clear nc, 
    orgsnclusters = onc.orgsnclusters;  clear onc,  
    trace_cluster = tc.trace_cluster; clear tc, 
    trace_x = tx.trace_x; clear tx;
    trace_y = ty.trace_y; clear ty;
    
    finite = get_finite_clusters(population,num_clusters,trace_cluster,trace_x,trace_y,INFRAT);
    clear trace_x trace_y trace_cluster,  
    %initialize correlation length holder
    monc = max(orgsnclusters(finite==1)); %maximum cluster mass
%     csd = zeros(monc,1);
    ngen = length(num_clusters);
    if ngen>TRANSIENCE,  start = TRANSIENCE+1;
    else, start = ngen;  end
    ns = zeros(numel(start:ngen),monc-limit+1); %generations x masses
    for gen = start:ngen, 
%       script_gen_update(this_script,gen,base_name,run_name);
      [cu,cv] = gen_index(num_clusters,gen); %generation indices
      onc = orgsnclusters(cu:cv); %masses in this generation
      [n,s] = hist(onc,limit:monc); %counted cluster sizes in range limit to max mass
%       norm = ceil(population(gen)/limit);
      ns(gen,1:(monc-limit+1)) = n./population(gen);%norm;
    end %for gen
    NS = cat_row(NS,ns);
%     [ns,s] = hist(orgsnclusters(cu:end),limit:monc);
%     NS = cat_row(NS,ns);
    end %trace_y
    end %trace_x
    end %trace_cluster
    end %orgsnclusters
    end %num_clusters
    end %population
  end %run
  cluster_number_distribution = mean(NS);
  % save in dir_name (\Data\) instead of new_dir_name (\Data\#_#\)
  save([dir_name 'cluster_number_distribution_' ...
    clus_name int2str(SIMS(1)) '_' int2str(SIMS(end))], ...
    'cluster_number_distribution');
  end %exists cluster_number_distribution
end %mu
end %dm
end %op      