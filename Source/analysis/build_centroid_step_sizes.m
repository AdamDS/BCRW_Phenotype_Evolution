%% build_centroid_step_sizes
%

this_script = 'build_centroid_step_sizes';
fprintf([this_script '\n']);
global SIMOPTS;

for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  [base_name,dir_name] = NameAndCD(make_dir,do_cd);
  clus_name = cluster_name(base_name);
  
  b = b +1;
  for run = SIMS, 
    run_name = int2str(run);
    new_dir_name = split_cd(dir_name,run,split,make_dir,do_cd);
    %% Load data
    [nc,go] = try_catch_load([new_dir_name 'num_clusters_' clus_name run_name],1,1);
    if go,  [cx,go] = try_catch_load([new_dir_name 'centroid_x_' clus_name run_name],1,1);
    if go,  [cy,go] = try_catch_load([new_dir_name 'centroid_y_' clus_name run_name],1,1);
    if go,  [ndc,go] = try_catch_load([new_dir_name 'num_descendant_clusters_' clus_name run_name],1,1);
    if go,  [dc,go] = try_catch_load([new_dir_name 'descendant_clusters_' clus_name run_name],1,1);
    if go,  
    num_clusters = nc.num_clusters; clear nc, 
    centroid_x = cx.centroid_x; clear cx, 
    centroid_y = cy.centroid_y; clear cy, 
    num_descendant_clusters = ndc.num_descendant_clusters;  clear ndc,  
    descendant_clusters = dc.descendant_clusters; clear dc,  
    
    %% Main algorithm
    ngen = length(num_clusters);
    centroid_step_sizes = size(descendant_clusters);
    for gen = 2:ngen, 
      [cu,cv] = gen_index(num_clusters,gen-1);  %parent cluster generation indices
      clusters = cu:cv; %parent cluster IDs
      cents = [centroid_x(clusters) centroid_y(clusters)];  %centroids of parent clusters
      [ccu,ccv] = gen_index(num_clusters,gen);  %child cluster generation indices
      child_clusters = ccu:ccv;  %child clusters IDs
      cc_cents = [centroid_x(child_clusters) centroid_y(child_clusters)]; %centroids of child clusters
      for c = clusters, %for each cluster
        num_dc = num_descendant_clusters(c);  %number of child clusters of c
        [dcu,dcv] = gen_index(num_descendant_clusters,c); %child cluster descendant indices
        desc_clus = descendant_clusters(dcu:dcv); %descendant clusters
        dc_cents = cc_cents(desc_clus,:); %centroids of desc_clus
        par_cent = cents(c,:);  %centroid of parent cluster
        dcents = sqrt(abs((par_cent -dc_cents).^2));  %calculate centroid step sizes
        centroid_step_sizes(desc_clus) = dcents; %store centroid step sizes
      end %for c
    end %for gen
    end %if descendant_clusters
    end %if num_descendant_clusters
    end %if centroid_y
    end %if centroid_x
    end %if num_clusters
    save([new_dir_name 'centroid_step_sizes_' clus_name run_name],'centroid_step_sizes');
  end %for run
end %for mu
end %for dm
end %for op