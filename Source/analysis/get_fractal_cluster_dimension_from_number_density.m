%% get_fractal_cluster_dimension_from_number_density.m
% Based on Feder's Fractals book p 34. Instead of particle density, this simply relates the number
% of organisms in a cluster
% ADS 3*8*13
global SIMOPTS;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  for run = SIMS
    run_name = int2str(run);
    go = 1;
    [p,go,error] = try_catch_load(['population_' base_name run_name],go,1);
    if go==1, [nc,go,error] = try_catch_load(['num_clusters_' base_name run_name],go,1);
    if go==1, [tx,go,error] = try_catch_load(['trace_x_' base_name run_name],go,1);
    if go==1, [ty,go,error] = try_catch_load(['trace_y_' base_name run_name],go,1);
    if go==1, [tc,go,error] = try_catch_load(['trace_cluster_' base_name run_name],go,1);
    if go==1, [tcs,go,error] = try_catch_load(['trace_cluster_seed_' base_name run_name],go,1);

    %initialize correlation length holder
    density = zeros(sum(num_clusters),1); %indexed by j
    j = 0;

    ngen = length(population);
    for gen = 1:ngen, 
      v = sum(population(1:gen)); u = v -population(gen) +1; 
      %get the cluster assignments of organisms in this generation
      tc_of_gen = trace_cluster(u:v);
      %get the coordinates of organisms in this generation
      coords_of_gen = [trace_x(u:v) trace_y(u:v)];
      %get nearest & second  nearest neighbor listings in this generation
      seeds_of_gen = trace_cluster_seed(u:v,:);
      
      %for each cluster
      for cluster = 1:num_clusters(gen)
        %get organisms of the same cluster
        these_orgs = find(tc_of_gen==cluster);
        %get x & y coords of only those in the same cluster
        these_coords = coords_of_gen(these_orgs,:);
        %get neighbors in the same cluster
        these_neighbors = neighbors_of_gen(these_orgs);
        %get alternates in the same cluster
        these_alternates = alternates_of_gen(these_orgs);
        distances_squared = zeros(sum(length(these_orgs)-1),1);
        i = 0;
        %for each organism except the last
        for org1 = 1:length(these_orgs(1:end-1))
          for org2 = 2:length(these_orgs(2:end))
            i = i +1;
            x1 = these_coords(org1,1);  y1 = these_coords(org1,2);
            x2 = these_coords(org2,1);  y2 = these_coords(org2,2);
            distances_squared(i) =  (x1-x2)^2 + (y1-y2)^2;
          end %org2
        end %org1
        j = j +1;
        correlation_lengths(j) = sqrt(max(distances_squared));
      end %cluster
    end %gen
  end %run
end %mu
end %dm
end %op      