%% get_infinite_probability.m *********************************
%
% -ADS 9*24*13
global SIMOPTS;
this_script = 'get_infinite_probability';
limit = SIMOPTS.limit;
tic;
%%
N = numel(overpop)*numel(death_max)*numel(mutability);
P = zeros(N,length(SIMS));
AVG_P = zeros(N,1);
STD_P = zeros(N,1);
i = 0;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  
  [base_name,dir_name] = NameAndCD(make_dir,do_cd);
  [clus_name] = cluster_name(base_name);
  
  fprintf('Attempting %s for %s \n',this_script,clus_name(1:end-1));
  
  i = i +1;
  for run = SIMS, 
    run_name = int2str(run);
    new_dir_name = split_cd(dir_name,run,split,make_dir,do_cd);

    r = find(run==SIMS);
    
    [p,go] = try_catch_load([new_dir_name 'population_' base_name run_name],1,1);
    if go,  [nc,go] = try_catch_load([new_dir_name 'num_clusters_' clus_name run_name],1,1);
    if go,  [tc,go] = try_catch_load([new_dir_name 'trace_cluster_' clus_name run_name],1,1);
    if go,  [tx,go] = try_catch_load([new_dir_name 'trace_x_' base_name run_name],1,1);
    if go,  [ty,go] = try_catch_load([new_dir_name 'trace_y_' base_name run_name],1,1);
    if go,  
      population = p.population;  clear p
      num_clusters = nc.num_clusters; clear nc
      trace_cluster = tc.trace_cluster; clear tc
      trace_x = tx.trace_x; clear tx
      trace_y = ty.trace_y; clear ty
      
      ngen = length(population);
      infinite = zeros(ngen,1);
%       if ngen>transience, start = transience; else, start = ngen; end
      try, INFRAT;  catch,  INFRAT = 0.989; end
      finite = get_finite_clusters(population,num_clusters,trace_cluster,...
        trace_x,trace_y,INFRAT);
%       for gen = 1:ngen, 
%         [cu,cv] = gen_index(num_clusters,gen);
%         if sum(finite(cu:cv))==num_clusters(gen), 
%           infinite(gen) = 0;
%         else, 
%           infinite(gen) = length(find(finite(cu:cv)<1));
%         end
%       end
%       P(i,r) = length(find(infinite~=0))./ngen;
      P(i,r) = length(find(finite~=1))./ngen;
    end
    end
    end
    end
    end
  end %for SIMS
  use = find(P(i,:)~=0);
  AVG_P(i) = mean(P(i,use));
  STD_P(i) = std(P(i,use));
end %mu
end %dm
end %op
toc;

if do_correlation_lengths_plot,  
  figure(525);
  errorbar(mutability,AVG_P,STD_P);
  xlabel('\mu');  ylabel('infinite cluster probability');
  title(make_title_name(generalize_base_name(clus_name),[]));
end