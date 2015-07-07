%% get_spatial_correlation_lengths.m
%
% Builds the data vector, correlation_lengths, for a set of data. The
% vector is organized like a centroids data vector, in that the order of
% the elements are blocked off by generation and by cluster identity order.
% The correlation lengths are taken to be the maximum spatial extent of
% each cluster. It does NOT correspond to the diameter of the cluster
% graph, which would be the greatest distance between related organisms
% within the cluster. That distance is needed for percolation correlation
% lengths, which may be the superior measure. These spatial correlation
% lengths may hopefully provide an alternative estimate of the same sort of
% quantity.
% -ADS 9*19*12
global SIMOPTS;
this_script = 'get_spatial_correlation_lengths';
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  for run = SIMS
    run_name = int2str(run);
    
    fprintf(['Attempting ' this_script ' for ' base_name run_name '\n']);
    if ~mat_exist(['correlation_lengths_' base_name run_name]) || SIMOPTS.write_over,  
    [p,go] = try_catch_load(['population_' base_name run_name],1,1); 
    if go,  [nc,go] = try_catch_load(['num_clusters_' base_name run_name],go,1); 
    if go,  [tx,go] = try_catch_load(['trace_x_' base_name run_name],go,1); 
    if go,  [ty,go] = try_catch_load(['trace_y_' base_name run_name],go,1); 
    if go,  [tc,go] = try_catch_load(['trace_cluster_' base_name run_name],go,1); 
    population = p.population;  clear p,  
    num_clusters = nc.num_clusters; clear nc, 
    trace_x = tx.trace_x; clear tx, 
    trace_y = ty.trace_y; clear ty, 
    trace_cluster = tc.trace_cluster; clear tc, 
    
    %initialize correlation length holder
    correlation_lengths = zeros(sum(num_clusters),1); %indexed by j
    j = 0;

    ngen = length(population);
    u = 0;  v = 0;
    for gen = 1:ngen
      script_gen_update(this_script,gen,base_name,run_name);
      u = v +1; v = sum(population(1:gen));
      %get the cluster assignments of organisms in this generation
      tc_of_gen = trace_cluster(u:v);
      %get the coordinates of organisms in this generation
      coords_of_gen = [trace_x(u:v) trace_y(u:v)];
      %for each cluster
      for cluster = 1:num_clusters(gen)
        %get organisms of the same cluster
        these_orgs = find(tc_of_gen==cluster);
        %get x & y coords of only those in the same cluster
        these_coords = coords_of_gen(these_orgs,:);
        distances_squared = zeros(sum(length(these_orgs)-1),1);
        i = 0;
        %for each organism except the last
        for org1 = 1:length(these_orgs(1:end-1))
          for org2 = 2:length(these_orgs(2:end))
            i = i +1;
            x1 = these_coords(org1,1);  y1 = these_coords(org1,2);
            x2 = these_coords(org2,1);  y2 = these_coords(org2,2);
            distances_squared(i) =  (x1-x2)^2 + (y1-y2)^2;
          end %for org2
        end %for org1
        j = j +1;
        correlation_lengths(j) = sqrt(max(distances_squared));
      end %for cluster
    end %for gen
    end %trace_cluster
    end %trace_y
    end %trace_x
    end %num_clusters
    end %exists correlation_lengths
  end %run
  save(['correlation_lengths_' base_name run_name],'correlation_lengths');
end %mu
end %dm
end %op      