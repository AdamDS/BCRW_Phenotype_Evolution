%% get_path_length.m
%This determines the distance of the path lentgh of clusters by summing seed_distances.
this_script = 'path_length';
fprintf([this_script '\n']);
global SIMOPTS;
% i = 0;
max_path_length_of_gen = zeros(NGEN*length(SIMS),total_params);
% last_max_path_length = zeros(num_final_gens,total_params);
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  r=0;
  for run = SIMS
    r = r+1;
    run_name = int2str(run);
    if exist([make_data_name('path_length',base_name,run_name,0)...
              '.mat'])~=2, 
    go = 1;
    [pop,go,error] = try_catch_load(['population_' base_name run_name],go,1);
    if go==1, [nc,go,error] = try_catch_load(['num_clusters_' base_name run_name],go,1);
    if go==1, [tc,go,error] = try_catch_load(['trace_cluster_' base_name run_name],go,1);
    if go==1, [sd,go] = try_catch_load(['seed_distances_' base_name run_name],go,1);
    if go==1,             
    fprintf([this_script ' for ' base_name run_name '\n']);
    
    population = pop.population;  clear pop
    num_clusters = nc.num_clusters; clear nc
    trace_cluster = tc.trace_cluster; clear tc
    seed_distances = sd.seed_distances; clear sd

    %initialize path_lengths
    path_length = zeros(sum(num_clusters),1);   
    
    ngen = size(num_clusters,2);
    k = 0;
    for gen = 1:ngen
      % the number of clusters in this generation
      script_gen_update(this_script,gen,base_name,run_name);
      
      v = sum(population(1:gen)); u = v -population(gen) +1; 
      cv = sum(num_clusters(1:gen));  cu = cv -num_clusters(gen) +1; 
      
      ref_pop = u:v;
      ref_clus = cu:cv;
      
      num_clus = num_clusters(gen);

      %initialize column vector for the per generation cluster path
      sd = seed_distances(ref_pop,:);        
      for c = 1:num_clus  %per generation loop calculates path lengths of cluster in one generation
        k = k +1;
        these = find(trace_cluster(ref_pop) == c);%finds the members of the ith cluster

        %seed_distances of jth cluster
        these_sd = sd(these,:); 
        path_length(k) = sum(sum(these_sd));
      end
      max_path_length_of_gen(gen,r) = max(path_length(cu:cv));
    end
%     i = i +1;
%     last_max_path_length(:,i) = max_path_length_of_gen(end-num_final_gens:end);
    save(['path_length_' base_name run_name],'path_length');
    end %seed_distances
    end %trace_cluster
    end %num_clusters
    end %population
    end %path_length 
  end %for run
end %for mu
end %for dm
end %for op

