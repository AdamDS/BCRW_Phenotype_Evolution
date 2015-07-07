%% get_correlation_lengths_distribution.m *********************************
% see Feder's Fractals p122 Eq7.22
% see Stauffer & Aharony's Introducion to Percolation Theory p60 Eq47b
% The gyration radii are based on distance between each organism and their
% cluster centroid, so the correlation lengths use Rg instead of 2Rg (this
% is used if the gyration radii are based on distance between organisms in
% the cluster).
% -ADS 5*1*13 updated 9*10*13
global SIMOPTS;
this_script = 'get_correlation_lengths_distribution';
limit = SIMOPTS.limit;
tic;
%%
N = numel(overpop)*numel(death_max)*numel(mutability);
avg_cl = zeros(N,length(SIMS));
AVG_CL = zeros(N,1);
STD_CL = zeros(N,1);

void = zeros(N,length(SIMS));
AVG_V = zeros(N,1);
STD_V = zeros(N,1);
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

    [cl,go] = try_catch_load([new_dir_name 'correlation_lengths_' clus_name run_name],1,1);
    if go,
      correlation_lengths = cl.correlation_lengths; clear cl
      ngen = length(correlation_lengths);
      if ngen>transience, start = transience; 
      else, start = floor(ngen/2);  end
      starts(i,r) = start;
      ngens(i,r) = ngen;
%       I = start:ngen;
%       use = isfinite(correlation_lengths(I)); %don't use NaN values
%       avg_cl(i,r) = mean(correlation_lengths(I(use==1)));
      avg_cl(i,r) = mean(correlation_lengths(start:ngen));
    end
    [p,gop] = try_catch_load([new_dir_name 'population_' base_name run_name],1,1);
    if gop,  
      population = p.population;  clear p
      A = 45*45;
      a = population*(op*op*pi/4);
      ngen = length(population);
      if ngen>transience, start = transience; 
      else, start = floor(ngen/2);  end
      void(i,r) = mean(A-a(start:ngen))/A;
    end
  end %for SIMS
  use = find(avg_cl(i,:)~=0);
  if ~isempty(use),
    defuse = find(avg_cl(i,use)~=NaN);
    AVG_CL(i) = mean(avg_cl(i,defuse));
    STD_CL(i) = std(avg_cl(i,defuse));
  end
  usep = find(void(i,:)~=0);
  if ~isempty(usep),  
    AVG_V(i) = mean(void(i,:));
    STD_V(i) = std(void(i,:));
  end
end %mu
end %dm
end %op
toc;

if do_correlation_lengths_plot,  
  figure(525);
  errorbar(mutability,AVG_CL,STD_CL);
  xlabel('\mu');  ylabel('\xi_{\perp}');
  title(make_title_name(generalize_base_name(clus_name),[]));
  
  figure(526);
  errorbar(mutability,AVG_V,STD_V);
  xlabel('\mu');  ylabel('\xi_{\perp}');
  title(['Void \xi_{\perp} ' make_title_name(generalize_base_name(clus_name),[])]);
end