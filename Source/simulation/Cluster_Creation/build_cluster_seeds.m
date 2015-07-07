%% build_cluster_seed - takes a string representing a simulation base
%function [trace_cluster_seed,seed_distances] = build_trace_cluster_seed(...
%  base_name,run,dir_name)
% The ith row of ouput corresponds to the ith organism as listed in
% trace_x and trace_y. The first column corresponds to the ith's nearest
% neighbor, and any subsequent column corresponds to the neighbors which
% are ranked in growing distance from i. For 3 limit runs (default), there
% are two columns, and for 2 limit, there is one column.
%
% IDEA: have a single output file which corresponds to the highest limit
% value run. For example, if 5 limit clusters are desired, then the
% information for 4 limit, 3 limit, and 2 limit, are still contained within
% the first however many rows. Only one file is needed in this case. There
% is a problem, however, if the simulations were run as 3 limit, and you
% want 3<limit clusters. If populations drop below the 3<limit, then a 4 or
% 5 limit cluster cannot be determined from the minimum population
% allowance of 3 limit.
%
% FIXED ON 22 AUG 2013; ALL trace_cluster_seed & seed_distances DATA 
% BEFORE for 2_limit MAY NEED TO BE FIXED!!!
% old: trace_cluster_seed = zeros(sum(population),2);
% new: trace_cluster_seed = zeros(sum(population),SIMOPTS.limit-1,'single');  %initialize tcs
%
% -ADS
function [trace_cluster_seed,seed_distances] = build_cluster_seeds(...
  base_name,run,dir_name)
global SIMOPTS;
this_function = 'build_cluster_seeds';

new_dir_name = split_cd(dir_name,run,SIMOPTS.split,1,0);
run_name = int2str(run);
clus_name = cluster_name(base_name);

print_function(this_function,[clus_name run_name]);

trace_cluster_seed = [];  seed_distances = [];
[clus_name] = cluster_name(base_name);
need_full = 1;
% if output files already exist for limit==3
if SIMOPTS.limit<3,
  if mat_exist([new_dir_name 'trace_cluster_seed_' base_name run_name]) && ...
     mat_exist([new_dir_name 'seed_distances_' base_name run_name]),  
   % load existing data, modify, and save
    [tcs,~,~] = try_catch_load([new_dir_name 'trace_cluster_seed_' base_name run_name]);
    [sd,~,~] = try_catch_load([new_dir_name 'seed_distances_' base_name run_name]);
    trace_cluster_seed = tcs.trace_cluster_seed;  clear tcs
    seed_distances = sd.seed_distances; clear sd
    trace_cluster_seed(:,2) = []; %shave off second nearest neighbors
    seed_distances(:,2) = []; %shave off second nearest neighbor distances
    save([new_dir_name 'trace_cluster_seed_' clus_name run_name],'trace_cluster_seed');
    save([new_dir_name 'seed_distances_' clus_name run_name],'seed_distances');
    need_full = 0;  %don't need to run the full function
  elseif mat_exist([new_dir_name 'trace_cluster_seed_' clus_name run_name]) || ...
     mat_exist([new_dir_name 'seed_distances_' clus_name run_name]), 
    need_full = 0;
  end
else, %limit==3
  if (mat_exist([new_dir_name 'trace_cluster_seed_' base_name run_name]) && ...
     mat_exist([new_dir_name 'seed_distances_' base_name run_name])) && ...
     ~SIMOPTS.write_over, %output files don't exist or write over
     need_full = 0;
  end %output files exist, so don't run
end

if need_full || SIMOPTS.write_over,
  [p,go] = try_catch_load([new_dir_name 'population_' base_name run_name],1,1); %get corresponding population

  if go,  
    if SIMOPTS.reproduction~=2, %if not random mating, need to measure distances
      [tx,go] = try_catch_load([new_dir_name ...
        'trace_x_' base_name run_name],1,1); %get corresponding trace_x
      [ty,go] = try_catch_load([new_dir_name ...
        'trace_y_' base_name run_name],1,1);  %get corresponding trace_x
    end
  if go,  
    population = p.population;  clear p
    trace_x = tx.trace_x; clear tx
    trace_y = ty.trace_y; clear ty

    NGEN = length(find(population)); %determine number of generations
    trace_cluster_seed = zeros(sum(population),SIMOPTS.limit-1,'single');  %initialize tcs
    seed_distances = zeros(sum(population),SIMOPTS.limit-1);
    
%FIXED ON 22 AUG 2013; ALL trace_cluster_seed & seed_distances DATA CREATED
%BEFORE NEEDS TO BE FIXED!!!
%     trace_cluster_seed = zeros(sum(population),2);

    ls = (((SIMOPTS.basic_map_size.*2)-1).*2)-1;
    warp_d = ls/2;
    u = 0; v = 0; %u and v index each population within trace_x/y
    for gen = 1:NGEN, 
      script_gen_update(this_function,gen,clus_name,run_name);
      u = v +1; %update lower limit index of this population
      v = sum(population(1:gen)); %update upper limit index of this population
      if SIMOPTS.reproduction~=2, 
        ix = trace_x(u:v);  %retrieve indiv x-coordinate for this population
        iy = trace_y(u:v);  %retrieve indiv y-coordinate for this population
        id = zeros(population(gen),1);
        for i = 1:population(gen), 
          if ~SIMOPTS.periodic(1) && ~SIMOPTS.periodic(2),  
            id = ((ix-ix(i)).^2) + ((iy-iy(i)).^2); %calculate distance from indiv i to all others
          else, %at least one periodic boundary
            if SIMOPTS.periodic(1), 
              adx = abs(ix(:)-ix(i)); %x distance from i to all others
              fixx = find(adx>warp_d(1)); %find those closer by wrap
              adx(fixx) = warp_d(1)-adx(fixx); %fix those closer by wrap
              dx2 = adx.^2;
            else, dx2 = (ix(:)-ix(i)).^2;  end
            if SIMOPTS.periodic(2), 
              ady = abs(iy(:)-iy(i)); %y distance from i to all others
              fixy = find(ady>warp_d(2)); %find those closer by wrap
              ady(fixy) = warp_d(2)-ady(fixy); %fix those closer by wrap
              dy2 = ady.^2;
            else, dy2 = (iy(:)-iy(i)).^2;  end
            id = dx2 +dy2;
          end
          sid = sort(id);
          %u+i-1 = index of current indiv i, adding u accounts for placement in
          %tcs while subtracting 1 is needed since u and i both start from 1 so
          %this adjustment is needed to obtain proper placement in tcs
          trace_cluster_seed(u+i-1,1) = find(id==sid(2)); %set mate to tcs
          seed_distances(u+i-1,1) = sid(2);
          if SIMOPTS.limit>=3, 
            trace_cluster_seed(u+i-1,2) = find(id==sid(3)); %set 2nd nearest neighbor to tcs
            seed_distances(u+i-1,2) = sid(3);
          end
        end
      elseif SIMOPTS.reproduction==2, %Random_Mating
        % Random Mating cluster seeds MUST be saved during simulation runtime.
        fprintf('ERROR: Did not save trace_cluster_seed at runtime for %s \n',[base_name run_name]);
      end %random mating
    end %gen
    seed_distances = sqrt(seed_distances);
    save([new_dir_name 'trace_cluster_seed_' clus_name run_name],'trace_cluster_seed');
    save([new_dir_name 'seed_distances_' clus_name run_name],'seed_distances');
  end %trace_x trace_y
  end %population
end %full function
end %function

%% Old Code
%         M = ceil(population(gen)*rand(population(gen),1));  %select a random mate for each organism
%         %check if mate and second nearest are the same, so need three unique indivs to make a 
%         %Nate cluster...maybe set SN post simulation......
%         same = find([1:population(gen)]'==M); %see if there are mates which identical to their partner
%         while length(same)>0, %while there is a mate that is itself
%           alt = ceil(rand(1,length(same))*population(gen)); %try an alternate mate
%           M(same) = alt; %set the alternative
%           same = find([1:population(gen)]'==M); %see if there are mates which are identical to their partner
%         end
%         SN = ceil(population(gen)*rand(population(gen),1));
%         same2 = find([1:population(gen)]'==SN | M==SN); %see if there are mates which identical to their partner
%         while length(same2)>0, %while there is a mate that is itself
%           alt = ceil(rand(1,length(same2))*population(gen)); %try an alternate mate
%           SN(same2) = alt; %set the alternative
%           same2 = find([1:population(gen)]'==SN | M==SN); %see if there are mates which are identical to their partner
%         end
%         trace_cluster_seed(u:v,:) = [M SN];
    %     clear M SN
