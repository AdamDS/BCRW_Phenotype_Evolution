%% get_corellation_length_growth_rate.m
%
i = 0;
N = numel(mutability)*numel(death_max)*numel(overpop);
gen_to_edge = inf(numel(SIMS),N);
tic;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  mineps = mu;
  maxeps = (2*(2*min(basic_map_size) -1) -1) -mu;
  i = i +1;
  for run = SIMS
    run_name = int2str(run);
    go = 1; [p,go] = try_catch_load(['population_' base_name run_name],go,1);
    if go==1, [x,go] = try_catch_load(['trace_x_' base_name run_name],go,1);
    if go==1, [y,go] = try_catch_load(['trace_y_' base_name run_name],go,1);
    if go==1, 
    population = p.population;  clear p
    trace_x = x.trace_x;  clear x
    trace_y = y.trace_y;  clear y
    at_edge = find(trace_x<mineps | trace_y<mineps | trace_x>maxeps | trace_y>maxeps);
    if at_edge, gen_to_edge(i) = at_edge(1);  end
    end %trace_y
    end %trace_x
    end %population
  end %for sims
end %for mu
end %for dm
end %for op