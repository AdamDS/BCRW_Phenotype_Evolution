%% build_correlation_lengths.m *********************************
% see Feder's Fractals p122 Eq7.22
% see Stauffer & Aharony's Introducion to Percolation Theory p60 Eq47b
% The gyration radii are based on distance between each organism and their
% cluster centroid, so the correlation lengths use Rg instead of 2Rg (this
% is used if the gyration radii are based on distance between organisms in
% the cluster).
%
% The default cluster size considered finite is less than 98.9% of a
% landscape side (INFRAT = 0.989).
% 
% The code within may be reused with other measures involving the cluster
% number distribution (within code, denoted as ns). I noted which lines (<)
% and blocks (v ... ^) are catered to the correlation lengths. 
%
% -ADS 5*1*13 updated 9*10*13
function [correlation_lengths] = build_correlation_lengths(base_name,run,dir_name), 
global SIMOPTS;
this_function = 'build_correlation_lengths'; % <
limit = SIMOPTS.limit;
INFRAT = 0.989;
correlation_lengths = [];

new_dir_name = split_cd(dir_name,run,SIMOPTS.split,1,0);
run_name = int2str(run);
[clus_name] = cluster_name(base_name);

print_function(this_function,[clus_name run_name]);

if ~mat_exist([new_dir_name 'correlation_lengths_' clus_name run_name]) || ... % <
    SIMOPTS.write_over,  
[p,go] = try_catch_load([new_dir_name 'population_' base_name run_name],1,1);
if go,  [tx,go] = try_catch_load([new_dir_name 'trace_x_' base_name run_name],1,1);
if go,  [ty,go] = try_catch_load([new_dir_name 'trace_y_' base_name run_name],1,1);
if go,  [nc,go] = try_catch_load([new_dir_name 'num_clusters_' clus_name run_name],1,1);
if go,  [tc,go] = try_catch_load([new_dir_name 'trace_cluster_' clus_name run_name],1,1);
if go,  [onc,go] = try_catch_load([new_dir_name 'orgsnclusters_' clus_name run_name],1,1);
if go,  [gr,go] = try_catch_load([new_dir_name 'gyration_radii_' clus_name run_name],1,1); % <
if go,  
  population = p.population;  clear p
  trace_x = tx.trace_x; clear tx
  trace_y = ty.trace_y; clear ty
  num_clusters = nc.num_clusters; clear nc
  trace_cluster = tc.trace_cluster; clear tc
  orgsnclusters = onc.orgsnclusters;  clear onc
  gyration_radii = gr.gyration_radii; clear gr % <

%   finite = get_finite_clusters(population,num_clusters,trace_cluster,trace_x,trace_y,INFRAT);
  clear trace_x trace_y trace_cluster

  ngen = length(find(population>=limit));
  monc = max(orgsnclusters);
  range = limit:monc;
  correlation_lengths = zeros(ngen,1);

  for gen = 1:ngen, 
    [cu,cv] = gen_index(num_clusters,gen);
%     f = finite(cu:cv); %f==1 is finite, f==0 is infinite along 1 axis, f==-1 is infinite both axes
    orgs = orgsnclusters(cu:cv);
    ORGS = [];
    if length(orgs)>1,  
      [~,toss] = max(orgs);
      I = 1:length(orgs); I(toss) = [];
      uo = unique(orgs(I));
      ORGS = orgs(I);
    else, 
      uo = unique(orgs);
      ORGS = orgs;  
    end
    n = hist(ORGS,range); %frequencies of finite s-clusters
    ns = n./population(gen); %cluster number distribution
    ssns = range(n~=0).*range(n~=0).*ns(n~=0); %ns*s^2
% v
    gr2 = gyration_radii(cu:cv).^2;
    avg_gr2 = zeros(size(uo))';
    j = 0;
    for s = uo, 
      j = j +1;
      if ~isempty(gr2(ORGS==s)),  
        avg_gr2(j) = mean(gr2(ORGS==s)); %<Rg^2(s)> of finite clusters
      end
    end
    correlation_lengths(gen) = sqrt((ssns*avg_gr2)./sum(ssns)); %sum(<Rg^2(s)>*ns*s^2)/sum(ns*s^2)
    clear avg_gr2
% ^
  end %for gen
  save([new_dir_name 'correlation_lengths_' clus_name run_name],'correlation_lengths'); % <
end %gyration_radii % <
end %orgsnclusters
end %trace_cluster
end %num_clusters
end %trace_y
end %trace_x
end %population
end %exist correlation_lengths % <
end %function