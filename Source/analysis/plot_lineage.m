% clear all;  clc;  close all;
%% plot_lineage.m
%
global SIMOPTS;
limit = SIMOPTS.limit;
n = 0;
N = numel(mutability)*numel(death_max)*numel(overpop);
tic;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  n = n +1;
  r = 0;
  for run = SIMS, 
    r = r +1;
    run_name = int2str(run);
    [p,go] = try_catch_load(['population_' base_name run_name],1,1);
    if go,  [nc,go] = try_catch_load(['num_clusters_' base_name run_name],go,1);
    if go,  [ncp,go] = try_catch_load(['num_clusters_produced_' base_name run_name],go,1);
    if go,  [cp,go] = try_catch_load(['clusters_produced_' base_name run_name],go,1);
    if go,  
    population = p.population;  clear p,  
    num_clusters = nc.num_clusters; clear nc, 
    num_clusters_produced = ncp.num_clusters_produced;  clear ncp,  
    clusters_produced = cp.clusters_produced; clear cp, 
    
    ngen = length(population);
    
    num_lins = num_clusters;
    num_pro = num_clusters_produced;
    pro = clusters_produced;
    M = max(num_clusters);
%     I = ngen +1;
    
%% Start debug level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear;  clc;  close all;
% ngen = 5;  %vertical time
% M = 6;  %horizontal lineage
% num_lins = [4 5 2 3 6]; %populations or clusters
% num_pro = [1 2 0 3, ...
%            2 0 0 1 0, ...
%            0 3, ...
%            3 1 2];  %number of forward time divergences
% pro = [1 2 3 3 4 5, ...
%        1 2 1, ...
%        1 2 3, ...
%        1 2 5 4 3 6];  %the forward time produce
% preN = [1 2 3 4, ...
%         1 2 3, ...
%         2 3 4, ...
%         2 3 4 5 6 7];
% connections = zeros(sum(num_pro),2);  %size of total produce x 2 (parent left, child right)
% connections = [1 1; ...
%                2 2; 2 3; ...
%                4 3; 4 4; 4 5; ...
%              1 1; 1 2; ...
%              3 1; ...
%            2 1; 2 2; 2 3; ...
%          1 1; 1 2; 1 5; ...
%          2 4; ...
%          3 3; 3 6]; 
% I = N +1;
%% Start main algorithm
%     scaled = [];
    a = 0;
    for i = 1:(length(num_pro)), 
      if num_pro(a),  
        for j = 1:num_pro(a), 
          a = a +1; 
          preN(a) = i
%       scaled = cat_row(scaled,[1:num_lins(i)].*M./(num_lins(i)));
    end
    figure(1), hold on;
    for i = 1:(ngen-1),
    %   I = I -1;
      nl = num_lins(i); %# lineages in ith generation               %4                        %5            
      v = sum(num_lins(1:i)); u = v -nl +1;                         %4, 1                     %9, 5         
      np = num_pro(u:v);  %# lineages produced                      %[1 2 0 3]                %[2 0 0 1 0]  
      vv = sum(num_pro(1:v)); uu = vv -sum(np) +1;                  %6, (6-6+1=)1             %9, (9-3+1=)7
      p = pro(uu:vv); %the lineages produced                        %[1 2 3 3 4 5]            %[1 2 1]      
    %   c = connections(uu:vv,:); %genealogy connections
      for n = 1:nl, %for each parent lineage
        if np(n), 
          vvv = sum(np(1:n)); uuu = vvv -np(n) +1;                  %1,1  %3,2    % %6,4      %8,7    % % %9,9
          par_of_n = uuu:vvv;%find(c(:,1)==n);                      %[1]  %[2 3]  % %[4 5 6]  %[7 8]  % % %[9]
    %     if par_of_n && par_of_n(1)~=,  %if any parent connections
    %       for pon = par_of_n(1):par_of_n(end), %for each parent connection
          chi_of_par_of_n = p(par_of_n);                            %[1]  %[2 3]  % %[3 4 5]  %[1 2]  % % %[1]
    %         chi_of_par_of_n = c(par_of_n,2);  %the child connections from parents
            for copon = chi_of_par_of_n(1):chi_of_par_of_n(end),  %for each child connection
%               plot([n copon],[i (i+1)]);%[I (I-1)]); 
%               plot([scaled(i,n) scaled(i+1,copon)],[i (i+1)]);%[I (I-1)]); 
              
              plot([n+preN copon+postN],[i (i+1)]);%[I (I-1)]); 
            end
    %       end
        else, 
          plot(scaled(i,n),i,'*k');%I,'*k');
        end
      end
    end
    hold off;
    xlim([0 M+1]);  xlabel('cluster lineage');  ylabel('generation');
%% End debug level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    title(make_title_name(base_name,run_name));
    end %clusters_produced
    end %num_clusters_produced
    end %num_clusters
    end %population
  end %for SIMS
end %for mu
end %for dm
end %for op