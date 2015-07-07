%% open_all_data.m
% Opens all babies data if it's available.
%
% -ADS 8*6*13
%

clus_name = cluster_name(base_name);

[p,gop,perror] = exist_load([new_dir_name 'population_' base_name run_name]);%,1,0);
[tx,gotx,txerror] = exist_load([new_dir_name 'trace_x_' base_name run_name]);%gop,0);
[ty,goty,tyerror] = exist_load([new_dir_name 'trace_y_' base_name run_name]);%gop,0);
[tcs,gotcs,tcserror] = exist_load([new_dir_name 'trace_cluster_seed_' clus_name run_name]);%gop,0);
[sd,gosd,sderror] = exist_load([new_dir_name 'seed_distances_' clus_name run_name]);%gotcs,0);
[par,gopar,parerror] = exist_load([new_dir_name 'parents_' base_name run_name]);%gop,0);
[k,gok,kerror] = exist_load([new_dir_name 'kills_' base_name run_name]);%gop,0);
[r,gor,rerror] = exist_load([new_dir_name 'rivalries_' base_name run_name]);%gop,0);
[nc,gonc,ncerror] = exist_load([new_dir_name 'num_clusters_' clus_name run_name]);%gop,0);
[onc,goonc,oncerror] = exist_load([new_dir_name 'orgsnclusters_' clus_name run_name]);%gonc,0);
[tc,gotc,tcerror] = exist_load([new_dir_name 'trace_cluster_' clus_name run_name]);%gonc,0);
[cx,gocx,cxerror] = exist_load([new_dir_name 'centroid_x_' clus_name run_name]);%gonc,0);
[cy,gocy,cyerror] = exist_load([new_dir_name 'centroid_y_' clus_name run_name]);%gonc,0);
[div,godiv,diverror] = exist_load([new_dir_name 'cluster_diversity_' clus_name run_name]);%gocx,0);
[dia,godia,diaerror] = exist_load([new_dir_name 'cluster_diameters_' clus_name run_name]);%gocx,0);
[gr,gogr,grerror] = exist_load([new_dir_name 'gyration_radii_' clus_name run_name]);%gocx,0);
[pl,gopl,plerror] = exist_load([new_dir_name 'path_length_' clus_name run_name]);%gocx,0);
[nd,gond,nderror] = exist_load([new_dir_name 'num_descendants_' base_name run_name]);%gop,0);
[ndc,gondc,ndcerror] = exist_load([new_dir_name 'num_descendant_clusters_' clus_name run_name]);%gotc,0);
[dc,godc,dcerror] = exist_load([new_dir_name 'descendant_clusters_' clus_name run_name]);%gondc,0);
[ncp,goncp,ncperror] = exist_load([new_dir_name 'num_clusters_produced_' clus_name run_name]);%gotc,0);
[cp,gocp,cperror] = exist_load([new_dir_name 'clusters_produced_' clus_name run_name]);%goncp,0);
[ncf,goncf,ncferror] = exist_load([new_dir_name 'num_clusters_fused_' clus_name run_name]);%gotc,0);
[cf,gocf,cferror] = exist_load([new_dir_name 'clusters_fused_' clus_name run_name]);%goncf,0);


% [p,gop,perror] = try_catch_load(['population_' base_name run_name],1,0);
% [tx,gotx,txerror] = try_catch_load(['trace_x_' base_name run_name],gop,0);
% [ty,goty,tyerror] = try_catch_load(['trace_y_' base_name run_name],gop,0);
% [tcs,gotcs,tcserror] = try_catch_load(['trace_cluster_seed_' clus_name run_name],gop,0);
% [sd,gosd,sderror] = try_catch_load(['seed_distances_' clus_name run_name],gotcs,0);
% [par,gopar,parerror] = try_catch_load(['parents_' base_name run_name],gop,0);
% [k,gok,kerror] = try_catch_load(['kills_' base_name run_name],gop,0);
% [r,gor,rerror] = try_catch_load(['rivalries_' base_name run_name],gop,0);
% [nc,gonc,ncerror] = try_catch_load(['num_clusters_' clus_name run_name],gop,0);
% [onc,goonc,oncerror] = try_catch_load(['orgsncluster_' clus_name run_name],gonc,0);
% [tc,gotc,tcerror] = try_catch_load(['trace_cluster_' clus_name run_name],gonc,0);
% [cx,gocx,cxerror] = try_catch_load(['centroid_x_' clus_name run_name],gonc,0);
% [cy,gocy,cyerror] = try_catch_load(['centroid_y_' clus_name run_name],gonc,0);
% [div,godiv,diverror] = try_catch_load(['cluster_diversity_' clus_name run_name],gocx,0);
% [dia,godia,diaerror] = try_catch_load(['cluster_diameters_' clus_name run_name],gocx,0);
% [pl,gopl,plerror] = try_catch_load(['path_length_' clus_name run_name],gocx,0);
% [nd,gond,nderror] = try_catch_load(['num_descendants_' base_name run_name],gop,0);
% [ndc,gondc,ndcerror] = try_catch_load(['num_descendant_clusters_' clus_name run_name],gotc,0);
% [dc,godc,dcerror] = try_catch_load(['descendant_clusters_' clus_name run_name],gondc,0);
% [ncp,goncp,ncperror] = try_catch_load(['num_clusters_produced_' clus_name run_name],gotc,0);
% [cp,gocp,cperror] = try_catch_load(['clusters_produced_' clus_name run_name],goncp,0);
% [ncf,goncf,ncferror] = try_catch_load(['num_clusters_fused_' clus_name run_name],gotc,0);
% [cf,gocf,cferror] = try_catch_load(['clusters_fused_' clus_name run_name],goncf,0);