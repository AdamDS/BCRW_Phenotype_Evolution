%% check_seed_distances.m
%
%
this_script = 'check_seed_distances';
fprintf([this_script '\n']);
global SIMOPTS;
N = length(overpop)*length(death_max)*length(mutability)*length(SIMOPTS.SIMS);
msd = zeros(N,SIMOPTS.limit-1);
totp = zeros(N,1);
bad_xy = [];
bad_sim = cell(N,1);
i = 0;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  for run = SIMS, 
    run_name = int2str(run);
    fprintf([this_script ' for ' base_name run_name '\n']);
    i = i +1;
    go = 1;
    [sd,go,error] = try_catch_load(['seed_distances_' base_name run_name],go,1);
    if go==1, [tcs,go,error] = try_catch_load(['trace_cluster_seed_' base_name run_name],go,1);
    if go==1, [p,go,error] = try_catch_load(['population_' base_name run_name],go,1);
    if go==1, [tx,go,error] = try_catch_load(['trace_x_' base_name run_name],go,1);
    if go==1, [ty,go,error] = try_catch_load(['trace_y_' base_name run_name],go,1);
    if go==1, 
      seed_distances = sd.seed_distances; clear sd
      trace_cluster_seed = tcs.trace_cluster_seed;  clear tcs
      population = p.population;  clear p
      trace_x = tx.trace_x; clear tx
      trace_y = ty.trace_y; clear ty
      
      if try_sqrt==1, seed_distances = sqrt(seed_distances);  end
      
      msd(i,:) = min(seed_distances((IPOP+1):end,:));
      bad_mates = find(seed_distances((IPOP+1):end,1)<op);
      bad_alts = find(seed_distances((IPOP+1):end,2)<op);
      totp(i) = sum(population);
      if min(msd(1,:))<op, 
        fprintf([make_data_name('seed_distances',base_name,run_name,0) ' is improper! %1.0f \n'],min(msd(1,:)));
      end
      if length(bad_mates)>0 || length(bad_alts)>0, 
        bad_xy = [bad_xy; [trace_x([bad_mates bad_alts]), trace_y([bad_mates bad_alts])]];
        bad_sim(i) = [base_name run_name];
      end
    end
    end
    end
    end
    end
  end
end
end
end