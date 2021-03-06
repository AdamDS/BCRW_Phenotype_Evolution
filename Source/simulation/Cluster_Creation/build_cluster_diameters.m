%% get_cluster_diameters.m
%
% Builds the data vector, cluster_diameters, for a set of data. The
% vector is organized like a centroids data vector, in that the order of
% the elements are blocked off by generation and by cluster identity order.
% -ADS 10*4*12
function [cluster_diameters] = build_diameters(base_name,run_name), 
this_script = 'get_cluster_diameters';
fprintf([this_script '\n']);
global SIMOPTS;
cluster_diameters = [];
[clus_name] = cluster_name(base_name);
if ~mat_exist(['cluster_diameters_' clus_name run_name]), 
[tcs,go,error] = try_catch_load(['trace_cluster_seed_' clus_name run_name],1,1);
if go, [tc,go,error] = try_catch_load(['trace_cluster_' clus_name run_name],go,1);
if go, [nc,go,error] = try_catch_load(['num_clusters_' clus_name run_name],go,1);
if go, [tx,go,error] = try_catch_load(['trace_x_' base_name run_name],go,1);
if go, [ty,go,error] = try_catch_load(['trace_y_' base_name run_name],go,1);
if go, [pop,go,error] = try_catch_load(['population_' base_name run_name],go,1);
if go,             
  fprintf(['get_cluster_diameters for ' clus_name run_name '\n']);
  trace_cluster_seed = tcs.trace_cluster_seed;  clear tcs
  trace_cluster = tc.trace_cluster;  clear tc
  num_clusters = nc.num_clusters;  clear nc
  trace_x = tx.trace_x; clear tx
  trace_y = ty.trace_y; clear ty
  population = pop.population;  clear pop

  limit = SIMOPTS.limit;
  %initialize characteristic length holders
  cluster_diameters = zeros(sum(num_clusters),1); %indexed by j
  ngen = length(population);
  u = 0;  v = 0;  cu = 0; cv = 0;
  for gen = 1:ngen, 
    script_gen_update(this_script,gen,clus_name,run_name);
    u = v +1; v = sum(population(1:gen));
    cu = cv +1; cv = sum(num_clusters(1:gen));
    %get the cluster assignments of organisms in this generation
    tc_of_gen = trace_cluster(u:v);
    %get the coordinates of organisms in this generation
    coords_of_gen = [trace_x(u:v) trace_y(u:v)];
    %get nearest & second  nearest neighbor listings in this generation
    relatives_of_gen = trace_cluster_seed(u:v,:);
    diameters = Inf*ones(num_clusters(gen),1);
    %for each cluster
    for cluster = 1:num_clusters(gen),  
%           if finite_clusters(cluster)==1
        %get organisms of the same cluster
        these_orgs = find(tc_of_gen==cluster);
        N = length(these_orgs);
        %get x & y coords of only those in the same cluster
        these_coords = coords_of_gen(these_orgs,:);
        %get seeds in the same cluster
        these_seeds = relatives_of_gen(these_orgs,:);
%% start of debugging
%         limit = 3; %2 if you want just nearest neighbors, 3 if you want nearest and alternate
%         these_orgs = [13; 15; 12; 11; 14; 16]; %sample list of organisms
%         N = length(these_orgs);
%         these_coords = [1 1; 1 3; 2 5; 9 9; 9 8; 9 1]; %coords for each organism
%         these_seeds = [15 12; 13 12; 15 13; 14 12; 11 12; 14 13]; %neighbors of these_orgs
%%
        %initialize shortest_paths matrix
        shortest_paths = Inf*ones(N,N);
        for i=1:N,  shortest_paths(i,i) = 0;  end
        %for each of these_orgs as the reference
        for reference_i = 1:length(these_orgs), 
          reference = these_orgs(reference_i);
          %initialize the list of organisms to_check
          to_check = these_orgs;
          checked = [];
          to_check_i = 1:length(these_orgs);
          %while there are still organisms to_check
          while ~isempty(to_check), 
            %get the shortest path to an organism not yet checked
            [from_d from] = min(shortest_paths(reference_i,to_check_i));
            %the organism that will be checked is from_this
            from_this = to_check(from); 
            from_i = find(these_orgs==from_this);
            %get coords of from_this
            x0 = these_coords(from_i,1);
            y0 = these_coords(from_i,2);
            %get connections
            to_them = these_seeds(from_i,:)'; %from_this to its seed
            to_them = unique([to_them; these_orgs(find(from_this==these_seeds(:,1)))]); %as a neighbor
            if limit>2, 
              to_them = unique([to_them; [to_them; these_orgs(find(from_this==these_seeds(:,2)))]]); %as an alternate
            end
            %remove any connections that have been checked
            remove = intersect(to_them,checked);
            for i = 1:length(remove)
              to_them(remove(i)==to_them) = [];
            end
            %if there is an organism to check in to_them
            if ~isempty(to_them), 
              %for each of the connections from_this to_them, check this_one
              for this_one = to_them',  
                this_i = find(this_one==these_orgs);
                %get coords of connections
                x1 = these_coords(this_i,1);
                y1 = these_coords(this_i,2);
                %measure distance from_this to its connections
                challenge_d_r_t = sqrt((x0-x1)^2 +(y0-y1)^2) +from_d;
                %get the old distances squared from reference to this_one
                old_d_r_t = shortest_paths(reference_i,this_i);
                %compare challenge distances to old distances
                [d_r_t picked_t] = min([challenge_d_r_t old_d_r_t]);
% Update shortest_paths considering latest results
                if picked_t==1, shortest_paths(reference_i,this_i) = d_r_t; end
              end
            end
% Update lists to_check & checked
            %remove this_one from to_check
            to_check(from_this==to_check) = []; checked = [checked; from_this];
            to_check_i(from_i==to_check_i) = [];
          end
        end
        diameters(cluster) = max(max(shortest_paths));
%% 
%           figure(914); plot(these_coords(:,1),these_coords(:,2),'x');
%           xlim([min(these_coords(:,1))-1 max(these_coords(:,1))+1]);
%           ylim([min(these_coords(:,2))-1 max(these_coords(:,2))+1]);
%% last part of debugging
    end %cluster
    cluster_diameters(cu:cv) = diameters;
  end %gen
  save(['cluster_diameters_' clus_name run_name],'cluster_diameters');
end
end
end
end
end
end
end
end

%% Old
% If a cluster is infinite (determined by INFRAT), then it has the value
% Inf. Therefore, all finite clusters are checked and infinite clusters are
% not considered. The default value for INFRAT is 0.9, which means that
% clusters which span at least 90% of the space in the x, y, or both axes,
% then it's cluster diameter is defaulted to Inf and it's cluster diameter
% is not determined by this script.

% for op = overpop, SIMOPTS.op = op;
% for dm = death_max, SIMOPTS.dm = dm;
% for mu = mutability, SIMOPTS.mu = mu;
%   make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
%   for run = SIMS
%     run_name = int2str(run);

  %       gyration_radii = zeros(sum(num_clusters),1);

  %         finite_clusters = get_finite_cluster_of_gen(tc_of_gen,coords_of_gen,INFRAT);

%               all_others(from_this==all_others) = [];

%           end %if finite_clusters

%       if INFRAT~=0.9
%         irn = num2str(INFRAT);
%         if INFRAT<1,  IRN = irn(3:end); else, IRN = [irn(1) '_' irn(3:end)];  end
%         save(['cluster_diameters_' IRN '_' base_name run_name],'cluster_diameters');
%       else
%       end

%   end %run
% end %mu
% end %dm
% end %op      