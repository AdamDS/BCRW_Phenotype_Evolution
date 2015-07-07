%% get_cluster_dimension.m
%
%
this_script = 'get_cluster_dimension';
fprintf([this_script '\n']);
global SIMOPTS;
init_z = zeros(length(overpop)*length(death_max)*length(mutability),length(SIMS));
M = init_z; B = init_z; EM = init_z;  EB = init_z;
SM = init_z; SB = init_z; SEM = init_z;  SEB = init_z;
max_cd = zeros(NGEN*length(SIMS),...
                              length(overpop)*length(death_max)*length(mutability));
i = 0;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  r = 0;  i = i +1;
  for run = SIMS, 
    j = 0;  r = r +1;
    run_name = int2str(run);
    fprintf([this_script ' for ' base_name run_name '\n']);
    go = 1;
    [p,go,error] = try_catch_load(['population_' base_name run_name],go,1);
    if go==1, [nc,go,error] = try_catch_load(['num_clusters_' base_name run_name],go,1);
    if go==1, [cdm,go,error] = try_catch_load(['cluster_diameters_' base_name run_name],go,1);
    if go==1, [onc,go,error] = try_catch_load(['orgsnclusters_' base_name run_name],go,1);
    if go==1, 
      population = p.population;  clear p
      num_clusters = nc.num_clusters; clear nc
      cluster_diameters = cdm.cluster_diameters;  clear cdm
      orgsnclusters = onc.orgsnclusters';  clear onc
      if length(cluster_diameters)>length(orgsnclusters), 
        fprintf('cd %1.0f > onc %1.0f \n',length(cluster_diameters),length(orgsnclusters));
      elseif length(cluster_diameters)<length(orgsnclusters), 
        fprintf('cd %1.0f < onc %1.0f \n',length(cluster_diameters),length(orgsnclusters));
      else, 
        fprintf('cd %1.0f = onc %1.0f \n',length(cluster_diameters),length(orgsnclusters));
      end
      ngen = length(population);
      if begin_of_end>ngen, start_gen = 4; 
      else, start_gen = begin_of_end; end
      for gen = start_gen:ngen, 
        j = j +1;
        v = sum(population(1:gen)); u = v -population(gen) +1;
        cv = sum(num_clusters(1:gen));  cu = cv -num_clusters(gen) +1;
        [max_cd(j,i),Imcd] = max(cluster_diameters(cu:cv));
        onc_gen = orgsnclusters(cu:cv);
        onc_of_max_cd(j,i) = onc_gen(Imcd);
      end %for gen
    end %orgsnclusters
    end %cluster_diameters
    end %num_clusters
    end %population
    r = r +1;
    M(i,r) = m;
    B(i,r) = b;
    EM(i,r) = em;
    EB(i,r) = eb;
    SM(i,r) = sm;
    SB(i,r) = sb;
    SEM(i,r) = sem;
    SEB(i,r) = seb;
  end %for run
end %for mu
end %for dm
end %for op