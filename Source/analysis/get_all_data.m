%% get_all_data.m
% get all babies data that is available
% 
% %Get data if available
% population = p.population;  clear p,  
% if gotx==1, trace_x = tx.trace_x; clear tx, end
% if goty==1, trace_y = ty.trace_y; clear ty,  end
% if gotcs==1,  trace_cluster_seed = tcs.trace_cluster_seed;  clear tcs,  end
% if gosd==1, seed_distances = sd.seed_distances; clear sd,  end
% if gopar==1,	parents = par.parents;  clear par,  end
% if gok==1,  kills = k.kills;  clear k,  end
% if gor==1,  rivalries = r.rivalries;  clear r,  end
% if gonc==1, num_clusters = nc.num_clusters; clear nc,  end
% if goonc==1,  orgsnclusters = onc.orgsnclusters;  clear onc,  end
% if gotc==1, trace_cluster = tc.trace_cluster; clear tc,  end
% if gocx==1, centroid_x = cx.centroid_x; clear cx,  end
% if gocy==1, centroid_y = cy.centroid_y; clear cy,  end
% if godiv==1,  cluster_diversity = div.cluster_diversity;  clear div,  end
% if godia==1,  cluster_diameters = dia.cluster_diameters;  clear dia,  end
% if gopl==1, path_length = pl.path_length; clear c, end
% if gond==1, num_descendants = nd.num_descendants; clear nd, end
% if gondc==1,  num_descendant_clusters = ndc.num_descendant_clusters;  clear ndc,  end
% if godc==1, descendant_clusters = dc.descendant_clusters; clear dc, end
% if goncp==1,  num_clusters_produced = ncp.num_clusters_produced;  clear ncp,  end
% if gocp==1, clusters_produced = cp.clusters_produced; clear cp, end
% if goncf==1,  num_clusters_fused = ncf.num_clusters_fused;  clear ncf,  end
% if gocf==1, clusters_fused = cf.clusters_fused; clear cf, end

%Get data if available
population = p.population;  clear p,  
if gotx, trace_x = tx.trace_x; clear tx, 
else, trace_x = []; end
if goty, trace_y = ty.trace_y; clear ty, 
else, trace_y = []; end
if gotcs,  trace_cluster_seed = tcs.trace_cluster_seed;  clear tcs,  
else, trace_cluster_seed = [];  end
if gosd, seed_distances = sd.seed_distances; clear sd,  
else, seed_distances = []; end
if gopar,	parents = par.parents;  clear par,  
else, parents = []; end
if gok,  kills = k.kills;  clear k,  
else, kills = []; end
if gor,  rivalries = r.rivalries;  clear r,  
else, rivalries = []; end
if gonc, num_clusters = nc.num_clusters; clear nc,  
else, num_clusters = []; end
if goonc,  orgsnclusters = onc.orgsnclusters;  clear onc,  
else, orgsnclusters = []; end
if gotc, trace_cluster = tc.trace_cluster; clear tc,  
else, trace_cluster = []; end
if gocx, centroid_x = cx.centroid_x; clear cx,  
else, centroid_x = []; end
if gocy, centroid_y = cy.centroid_y; clear cy,  
else, centroid_y = []; end
if godiv,  cluster_diversity = div.cluster_diversity;  clear div,  
else, cluster_diversity = []; end
if godia,  cluster_diameters = dia.cluster_diameters;  clear dia,  
else, cluster_diameters = []; end
if gogr,  gyration_radii = gr.gyration_radii;  clear gr,  
else, gyration_radii = []; end
if gopl, path_length = pl.path_length; clear c,  
else, path_length = []; end
if gond, num_descendants = nd.num_descendants; clear nd,  
else, num_descendants = []; end
if gondc,  num_descendant_clusters = ndc.num_descendant_clusters;  clear ndc,  
else, num_descendant_clusters = []; end
if godc, descendant_clusters = dc.descendant_clusters; clear dc,  
else, descendant_clusters = []; end
if goncp,  num_clusters_produced = ncp.num_clusters_produced;  clear ncp,  
else, num_clusters_produced = []; end
if gocp, clusters_produced = cp.clusters_produced; clear cp,  
else, clusters_produced = []; end
if goncf,  num_clusters_fused = ncf.num_clusters_fused;  clear ncf,  
else, num_clusters_fused = []; end
if gocf, clusters_fused = cf.clusters_fused; clear cf,  
else, clusters_fused = []; end