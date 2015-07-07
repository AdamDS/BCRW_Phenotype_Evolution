%% build_cluster_numbers.m
%
%
global SIMOPTS;
this_script = 'build_cluster_numbers';
limit = SIMOPTS.limit;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  [base_name,dir_name] = NameAndCD(0,0);
  [clus_name] = cluster_name(base_name);
  
  fprintf('Attempting %s for %s \n',this_script,clus_name(1:end-1));
  
  cn_name = make_cluster_numbers_name(dir_name,clus_name,SIMS);
  if ~mat_exist(cn_name) || SIMOPTS.write_over,  
    for run = SIMS, 
      run_name = int2str(run);

      new_dir_name = split_cd(dir_name,run,split,0,0);

      fprintf(' %s of %s \n',run_name,int2str(length(SIMS)));

      [mf,go] = try_catch_load([new_dir_name 'mass_frequencies_' clus_name run_name],1,1); 
      if go,  
        mass_frequencies = mf.mass_frequencies; clear mf
        MF = cat_row(MF,mass_frequencies);
      end %mass_frequencies
    end %run
    cluster_numbers = mean(MF,1);
  end %exists cluster_number_distribution
  save(cn_name,'cluster_numbers'); 
end %mu
end %dm
end %op      