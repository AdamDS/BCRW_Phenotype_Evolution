%% initialize_check_values.m
% initializes check logicals for all babies data
%
% %Checks on populaion
% c_ipop = 0; %initial populaion
% c_ngen = 0; %ngen
% c_minp = 0; %minimum population
% 
% %Checks on data size
% c_lentx = 0;  %rows of trace_x
% c_lenty = 0;  %rows of trace_y
% c_lentcs = 0; %rows of trace_cluster_seed
% c_lensd = 0;  %rows of seed_distances
% c_lenk = 0; %rows of kills
% c_lenr = 0; %rows of rivalries
% c_lenpar = 0; %rows of parents
% c_lennc = 0;  %length of num_clusters
% c_lentc = 0;  %rows of trace_cluster
% c_lenonc = 0; %rows of orgsncluster
% c_lencx = 0;  %rows of centroid_x
% c_lency = 0;  %rows of centroid_y
% c_lendiv = 0; %rows of cluster_diversity
% c_lendia = 0; %rows of cluster_diameters
% c_lenpl = 0;  %rows of path_length
% c_rownd = 0;  %rows of num_descendants
% c_colnd = 0;  %columns of num_descendants
% c_rowndc = 0; %rows of num_descendant_clusters
% c_colndc = 0; %columns of num_descendant_clusters
% c_lendc = 0;  %rows of descendant_clusters
% c_lenncp = 0; %rows of num_clusters_produced
% c_lencp = 0;  %rows of clusters_produced
% c_lenncf = 0; %rows of num_clusters_fused
% c_lencf = 0;  %rows of clusters_fused
% 
% %Checks on trace_x
% c_mintx = 0;  %minimum trace_x
% c_maxtx = 0;  %maximum trace_x
% 
% %Checks on trace_y
% c_minty = 0;  %minimum trace_y
% c_maxty = 0;  %maximum trace_y
% 
% %Checks on trace_cluster_seed
% c_mintcs = 0; %minimum trace_cluster_seed
% c_maxtcs = 0; %maximum trace_cluster_seed
% 
% %Checks on seed_distances
% c_minsd = 0;  %minimum seed_distances
% 
% %Checks on kills
% c_maxk = 0; %maximum kills
% c_maxsumk = 0;  %maximum sum of all kills
% 
% %Checks on rivalries
% c_maxr = 0; %maximum rivalries
% 
% %Checks on parents
% c_minpar = 0; %minimum parents
% c_maxpar = 0; %maximum parents
% 
% %Checks on num_clusters
% c_minnc = 0;  %minimum num_clusters
% c_maxnc = 0;  %maximum num_clusters
% 
% %Checks on trace_cluster
% c_mintc = 0;  %minimum trace_cluster
% c_maxtc = 0;  %maximum trace_cluster
% 
% %Checks on orgsncluster
% c_minonc = 0; %minimum orgsncluster
% c_maxonc = 0; %maximum orgsncluster
% 
% %Checks on centroid_x
% c_mincx = 0;  %minimum centroid_x
% c_maxcx = 0;  %maximum centroid_x
% 
% %Checks on centroid_y
% c_mincy = 0;  %minimum centroid_y
% c_maxcy = 0;  %maximum centroid_y
% 
% %Checks on cluster_diversity
% c_mindiv = 0; %minimum cluster_diversity
% c_maxdiv = 0; %maximum cluster_diversity
% 
% %Checks on cluster_diameters
% c_mindia = 0; %minimum cluster_diameters
% 
% %Checks on path_length
% c_minpl = 0;  %minimum path_length
% 
% %Checks on num_descendants
% c_minnd = 0;  %minimum num_descendants
% c_maxnd = 0;  %maximum num_descendants
% 
% %Checks on num_descendant_clusters
% c_minndc = 0; %minimum num_descendant_clusters
% c_maxndc = 0; %maximum num_descendant_clusters
% 
% %Checks on descendant_clusters
% c_mindc = 0;  %minimum descendant_clusters
% c_maxdc = 0;  %maximum descendant_clusters
% 
% %Checks on num_clusters_produced
% c_minncp = 0; %minimum num_clusters_produced
% c_maxncp = 0; %maximum num_clusters_produced
% 
% %Checks on clusters_produced
% c_mincp = 0;  %minimum clusters_produced
% c_maxcp = 0;  %maximum clusters_produced
% 
% %Checks on num_clusters_fused
% c_minncf = 0; %minimum num_clusters_fused
% c_maxncf = 0; %maximum num_clusters_fused
% 
% %Checks on clusters_fused
% c_mincf = 0;  %minimum clusters_fused
% c_maxcf = 0;  %maximum clusters_fused

%Checks on populaion
c_ipop = 0; %initial populaion
c_ngen = 0; %ngen
c_minp = 0; %minimum population

%Checks on data size
c_lentx = 0;  %rows of trace_x
c_lenty = 0;  %rows of trace_y
c_lentcs = 0; %rows of trace_cluster_seed
c_lensd = 0;  %rows of seed_distances
c_lenk = 0; %rows of kills
c_lenr = 0; %rows of rivalries
c_lenpar = 0; %rows of parents
c_lennc = 0;  %length of num_clusters
c_lentc = 0;  %rows of trace_cluster
c_lenonc = 0; %rows of orgsncluster
c_lencx = 0;  %rows of centroid_x
c_lency = 0;  %rows of centroid_y
c_lendiv = 0; %rows of cluster_diversity
c_lendia = 0; %rows of cluster_diameters
c_lenpl = 0;  %rows of path_length
c_rownd = 0;  %rows of num_descendants
c_colnd = 0;  %columns of num_descendants
c_rowndc = 0; %rows of num_descendant_clusters
c_colndc = 0; %columns of num_descendant_clusters
c_lendc = 0;  %rows of descendant_clusters
c_lenncp = 0; %rows of num_clusters_produced
c_lencp = 0;  %rows of clusters_produced
c_lenncf = 0; %rows of num_clusters_fused
c_lencf = 0;  %rows of clusters_fused

%Checks on trace_x
c_mintx = 0;  %minimum trace_x
c_maxtx = 0;  %maximum trace_x

%Checks on trace_y
c_minty = 0;  %minimum trace_y
c_maxty = 0;  %maximum trace_y

%Checks on trace_cluster_seed
c_mintcs = 0; %minimum trace_cluster_seed
c_maxtcs = 0; %maximum trace_cluster_seed

%Checks on seed_distances
c_minsd = 0;  %minimum seed_distances

%Checks on kills
c_maxk = 0; %maximum kills
c_maxsumk = 0;  %maximum sum of all kills

%Checks on rivalries
c_maxr = 0; %maximum rivalries

%Checks on parents
c_minpar = 0; %minimum parents
c_maxpar = 0; %maximum parents

%Checks on num_clusters
c_minnc = 0;  %minimum num_clusters
c_maxnc = 0;  %maximum num_clusters

%Checks on trace_cluster
c_mintc = 0;  %minimum trace_cluster
c_maxtc = 0;  %maximum trace_cluster

%Checks on orgsncluster
c_minonc = 0; %minimum orgsncluster
c_maxonc = 0; %maximum orgsncluster

%Checks on centroid_x
c_mincx = 0;  %minimum centroid_x
c_maxcx = 0;  %maximum centroid_x

%Checks on centroid_y
c_mincy = 0;  %minimum centroid_y
c_maxcy = 0;  %maximum centroid_y

%Checks on cluster_diversity
c_mindiv = 0; %minimum cluster_diversity
c_maxdiv = 0; %maximum cluster_diversity

%Checks on cluster_diameters
c_mindia = 0; %minimum cluster_diameters

%Checks on path_length
c_minpl = 0;  %minimum path_length

%Checks on num_descendants
c_minnd = 0;  %minimum num_descendants
c_maxnd = 0;  %maximum num_descendants

%Checks on num_descendant_clusters
c_minndc = 0; %minimum num_descendant_clusters
c_maxndc = 0; %maximum num_descendant_clusters

%Checks on descendant_clusters
c_mindc = 0;  %minimum descendant_clusters
c_maxdc = 0;  %maximum descendant_clusters

%Checks on num_clusters_produced
c_minncp = 0; %minimum num_clusters_produced
c_maxncp = 0; %maximum num_clusters_produced

%Checks on clusters_produced
c_mincp = 0;  %minimum clusters_produced
c_maxcp = 0;  %maximum clusters_produced

%Checks on num_clusters_fused
c_minncf = 0; %minimum num_clusters_fused
c_maxncf = 0; %maximum num_clusters_fused

%Checks on clusters_fused
c_mincf = 0;  %minimum clusters_fused
c_maxcf = 0;  %maximum clusters_fused