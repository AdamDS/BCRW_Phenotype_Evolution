%% get_same_cluster_correlation_length_distribution.m ***************************
global SIMOPTS;
this_script = 'get_same_cluster_correlation_length_distribution';
tic;
%defaults
bms = SIMOPTS.basic_map_size;
X = (((2*bms(1))-1)*2)-1; Y = (((2*bms(2))-1)*2)-1;
max_land = sqrt(X^2 +Y^2);
regions = 50;
coeps = max_land/regions;
cohoods = [0:coeps:max_land];
transient = 201;
%% initialize outputs
%%
for op = overpop, SIMOPTS.op = op;  ophoods = [op:op:max_land];
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  [base_name,~] = NameAndCD(0,do_cd);
  for run = SIMS
    run_name = int2str(run);
    fprintf([this_script ' for ' base_name run_name '\n']);
    [p,go] = try_catch_load(['population_' base_name run_name],1,1);
    if go,  [x,go] = try_catch_load(['trace_x_' base_name run_name],go,1);
    if go,  [y,go] = try_catch_load(['trace_y_' base_name run_name],go,1);
    if go,  [nc,go] = try_catch_load(['num_clusters_' base_name run_name],go,1);
    if go,  [tc,go] = try_catch_load(['trace_cluster_' base_name run_name],go,1);
    if go,
      population = p.population;  clear p
      trace_x = x.trace_x;  clear x
      trace_y = y.trace_y;  clear y
      num_clusters = nc.num_clusters; clear nc
      trace_cluster = tc.trace_cluster; clear tc
      ngen = length(population);
      for gen = transient:ngen,
        script_gen_update(this_script,gen,base_name,run_name);
        [u,v] = gen_index(population,gen);
        [cu,cv] = gen_index(num_clusters,gen);
        x = trace_x(u:v);
        y = trace_y(u:v);
        c = trace_cluster(cu:cv);
        all_dists = pdist([x y]);
        avg_dists(gen) = mean(all_dists);
      end %for gen
    end %trace_cluster
    end %num_clusters
    end %trace_y
    end %trace_x
    end %population
  end %for run
end %for mu
end %for dm
end %for op
toc;

if do_characteristic_diameters==1, 
  [asdf,qwer] = hist(CD_finite,100);
  [m,b,sm,sb] = linear_fit(log10(Neighborhoods(25:75)'),log10((25:75)),log10((25:75)),1);
end
if do_correlation_lengths==1, 
  [mcl,bcl,smcl,sbcl] = linear_fit(log10(Neighborhoods(5:15)'),log10(AVG_N(5:15,end)),log10(STD_N(5:15)),1);
end
% xlswrite(['characteristic_lengths_' generalize_base_name(base_name) '__5sim'],...
%          [Neighborhoods, AVG_N]);