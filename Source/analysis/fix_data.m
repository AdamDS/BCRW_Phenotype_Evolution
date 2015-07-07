%% fix_data.m
%
%
this_script = 'fix_data';
fprintf([this_script '\n']);
global SIMOPTS;
N = length(overpop)*length(death_max)*length(mutability)*length(SIMOPTS.SIMS);
limit = SIMOPTS.limit;
min_land = 0.5;
max_land = 0.5 +2*((2*max(basic_map_size)) -1) -1;
smax_land = 0.5 +2*((2*min(basic_map_size)) -1) -1;
diag_land = sqrt(max_land^2 +smax_land^2);
i = 0;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  for run = SIMS, 
    run_name = int2str(run);
    fprintf([this_script ' for ' base_name run_name '\n']);
    i = i +1;
    if do_fix_population, fix_population(base_name,run_name); end
    if do_fix_gyration_radii, fix_gyration_radii(base_name,run,dir_name); end
  end
end
end
end