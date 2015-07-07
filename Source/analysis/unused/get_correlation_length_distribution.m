%% get_correlation_length_distribution.m *********************************
%
% -ADS 5*1*13
global SIMOPTS;
tic;
%defaults
%% initialize outputs
%%
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir,do_cd);
%   i = find(mu==mutability);
  if ~mat_exist(['correlation_lengths_' base_name(1:end-1)]), 
  for run = SIMS, 
    run_name = int2str(run);
    [cnd,go] = try_catch_load(['cluster_number_distribution_' base_name run_name],1,1);
    if go,  [gr,go] = try_catch_load(['gyration_radii_' base_name run_name],go,1);
    if go,  [onc,go] = try_catch_load(['orgsnclusters_' base_name run_name],go,1);
    if go,  
    cluster_number_distribution = cnd.cluster_number_distribution;  clear cnd
    gyration_radii = gr.gyration_radii; clear gr
    orgsnclusters = onc.orgsnclusters;  clear onc
    for s=limit:length(cluster_number_distribution),  
      ns = cluster_number_distribution(s);
      is = find(orgsnclusters==s);
      Rg = mean(gyration_radii(is));
      num(s) = 2*Rg*(s^2)*ns;
      den(s) = (s^2)*ns;
    end % for s
    correlation_lengths = sqrt(sum(num)./sum(den));
    end % orgsnclusters
    end % gyration_radii
    end % cluster_number_distribution
  end % for SIMS
  if go,  save(['correlation_lengths_' base_name(1:end-1)],'correlation_lengths'); end
  end % exists correlation_lengths
end
end
end
toc;

