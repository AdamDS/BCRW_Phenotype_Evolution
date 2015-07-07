clear;  clc;
%% Debug level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_clusters = [2 6 3]; %indexed by generation
centroid_x = [3 7, 1 2 3 4 5 6, 2 5 8]';
centroid_y = [2 9, 2 6 1 4 8 5, 4 2 7]';
orgsnclusters = [3 5, 4 5 7 3 5 4, 4 3 3]';
num_clusters_produced = [3 5, 1 1 2 1 2 1]'; %indexed for num_clusters(2:end) 
                                              %[sum >= num_clusters(2), sum >= num_clusters(3), ...]
clusters_produced = [1 2 3,1 2 4 5 6, 1,1,1 2,1,2 3,1]'; %indexed for sum(num_clusters_produced)
dbg_tree_links = [1;1;1;2;2;2;2;2; 1;2;3;3;4;5;5;6, clusters_produced];
dbg_dcent_links = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
ngen = length(num_clusters);
cu = 0; cv = 0; %indexes current generation cluster data
cpu = 0;  cpv = 0;  %indexes clusters produced data
icu = 0;  icv = 0;  %indexes output
tree_links = zeros(length(clusters_produced),2); % [parent cluster, child cluster]
dcent = zeros(length(clusters_produced),1); % distance between parent cluster & child cluster
dcent_links = tree_links;
prop = zeros(length(clusters_produced),1);  % ratio of cluster size, child to parent cluster
prop_links = tree_links;
%% get the clusters of this gen
for gen = 1:ngen-1, 
  cu = cv +1; cv = sum(num_clusters(1:gen));  %indexes this gen stuff
  cents_of_gen = [centroid_x(cu:cv), centroid_y(cu:cv)];  %centroids of this gen
  orgs_of_gen = orgsnclusters(cu:cv); %orgs of this gen
  ncps_of_gen = num_clusters_produced(cu:cv); %children produced of this gen
  
  cpu = cpv +1; cpv = sum(num_clusters_produced(1:cv)); %indexes next gen produced stuff
  cp_next_gen = clusters_produced(cpu:cpv); %children clusters of next gen
  gbeg = cv +1;  gend = sum(num_clusters(1:gen+1)); %indexes next gen stuff
  cents_next_gen = [centroid_x(gbeg:gend), centroid_y(gbeg:gend)];  %centroids of next gen
  orgs_next_gen = orgsnclusters(gbeg:gend); %orgs of next gen
  
  cbeg = 0; cend = 0;
  % for each cluster
  for c = 1:num_clusters(gen), 
    cbeg = cend +1; cend = sum(ncps_of_gen(1:c));
    % get number of children clusters
    nprod = ncps_of_gen(c);
    icu = icv +1; icv = icu +nprod;
    % get the children clusters
    cp = cp_next_gen(cbeg:cend);
    tree_links(icu:icv,:) = [c*ones(nprod,1), cp];
    % get the centroid of the this cluster
    cent = cents_of_gen(c,:);
    % get the centroids of each child cluster
    ccent = cents_next_gen(cp,:);
    % sort children clusters based on centroid distance from parent cluster
    [sdbranches, idb] = sort((cent(1) -ccent(:,1)).^2 +(cent(2) -ccent(:,2)).^2); %s(i) = d
    % assign a location of the children clusters 
    dcent(icu:icv) = sdbranches;
    dcent_links(icu:icv,:) = idb;
    % get the number of organisms of each child cluster
    corg = orgs_next_gen(cp,:);
    % sort children clusters based on number of organisms they contain
    
  end
end