% function [nearest_neighbor_distances] = get_nearest_neighbor_distances(base_name,run_name,limit,do_parents_plot)
% exp_name = [base_name run_name];
% pn = ['population_' exp_name];  load(pn);
% tcsn = ['trace_cluster_seed_' exp_name];  load(tcsn);
% txn = ['trace_x_' exp_name];  load(txn);
% tyn = ['trace_y_' exp_name];  load(tyn);
% nearest_neighbor_distances = zeros(sum(population),1);
% u = 1;  v = 0;
% for gen = 1:length(population)
%   v = sum(population(1:gen));
%   if population(gen)>limit
%     refs = u:v;
%     mates = u +trace_cluster_seed(u:v,1) -1;
%     nearest_neighbor_distances(u:v,1) = sqrt((trace_x(refs)-trace_x(mates)).^2 + (trace_y(refs)-trace_y(mates)).^2);
%   end
%   u = v +1;
% end
% if do_parents_plot==1
%   hist(nearest_neighbor_distances,ceil(sum(population)*0.01));
% end
% end
this_script = 'get_R';
fprintf([this_script '\n']);
global SIMOPTS;
% max_size = 3;%(((max(basic_map_size)*2)-1)*2)-1;
AVG_GEN_NN_DIST = zeros(NGEN,length(mutability));
STD_GEN_NN_DIST = zeros(NGEN,length(mutability));
avg_nn_dist = zeros(length(mutability),length(SIMS));
AVG_NN_DIST = zeros(length(mutability),1);
STD_NN_DIST = zeros(length(mutability),1);

AVG_GEN_NC_DIST = zeros(NGEN,length(mutability));
STD_GEN_NC_DIST = zeros(NGEN,length(mutability));
avg_nc_dist = zeros(length(mutability),length(SIMS));
AVG_NC_DIST = zeros(length(mutability),1);
STD_NC_DIST = zeros(length(mutability),1);

i = 1;
avgR = zeros(length(mutability),length(SIMS));
avgc = zeros(length(mutability),length(SIMS));
AVG_R = zeros(length(mutability),1);
STD_R = zeros(length(mutability),1);
AVG_c = zeros(length(mutability),1);
STD_c = zeros(length(mutability),1);
GEN_R = zeros(NGEN,length(mutability));
GEN_c = zeros(NGEN,length(mutability));

b = 0;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  b = b +1;  
  m_nn_dist = zeros(NGEN,length(SIMS));
  r = 0;
  for run = SIMS, 
    r = r +1;
    run_name = int2str(run);
    new_dir_name = split_cd(dir_name,run,split,make_dir,do_cd);
    clus_name = cluster_name(base_name);
    long = 0;
    %need population, trace_cluster_seed, and seed_distances as minimum requirements
    [p,go] = try_catch_load([new_dir_name 'population_' base_name run_name],1,1);
    if go,  %have population
      population = p.population;  clear p,  
      [tcs,go] = try_catch_load([new_dir_name 'trace_cluster_seed_' clus_name run_name],1,1);
      if go,  trace_cluster_seed = tcs.trace_cluster_seed; clear tcs, end
    end
    if go,  %have trace_cluster_seed
      [sd,go] = try_catch_load([new_dir_name 'seed_distances_' clus_name run_name],1,1);
      if go,  seed_distances = sd.seed_distances; clear sd, end
    end
    
    if do_population_level, 
      if go, %have seed_distances
        fprintf([this_script ' for populations of ' base_name run_name '\n']);
        [tx,go,txe] = try_catch_load([new_dir_name 'trace_x_' base_name run_name],1,1);
        if go,  [ty,go,txe] = try_catch_load([new_dir_name 'trace_y_' base_name run_name],1,1);
        if go, 
          trace_x = tx.trace_x; clear tx
          trace_y = ty.trace_y; clear ty

          ngen = length(population);

          lsx = (((basic_map_size(2)*2) -1)*2) -1;
          lsy = (((basic_map_size(1)*2) -1)*2) -1;
          u = 0;  v = 0;
          for gen = 1:ngen %for each generation of this run
            u = v +1;
            v = sum(population(1:gen));
            refs = u:v; %the reference organisms of this generation
            m_nn_dist(gen,run) = mean(seed_distances(refs,1));  %average nearest neighbor distances
          end %for gen
          this_m_nn_dist = m_nn_dist(find(m_nn_dist(:,r)),r); %m_nn_dist of this simulation
          avg_nn_dist(b,r) = mean(this_m_nn_dist); %average of averages nn_dists
          %%
          area = lsx*lsy; %area of the landscape
          rho = population'./area;  %density for each generation
          mrd = 0.5*(rho.^-0.5);  %denominator of R (expected density if random population)
          R = this_m_nn_dist./mrd;  %R for each generation
          sigma_mrd = 0.26136./sqrt(population'.*rho);  %variance of R for each generation
          c = (this_m_nn_dist -mrd)./sigma_mrd; %confidence in R for each generation
          avgR(b,run) = mean(R);  %time averaged R
          avgc(b,run) = mean(c);  %time averaged c
          GEN_R(1:length(R),b) = R; %store R for each generation for this simulation run
          GEN_c(1:length(c),b) = c; %store c for each generation for this simulation run
    %       AVG_GEN_NN_DIST(:,b) = mean(this_m_nn_dist);  %
    %       STD_GEN_NN_DIST(:,b) = std(this_m_nn_dist);
        end
        end
      end
        %%
    elseif do_cluster_level==1, 
      fprintf([this_script ' for clusters of ' clus_name run_name '\n']);
      [nc,go,~] = try_catch_load([new_dir_name 'num_clusters_' clus_name run_name],1,1);
      if go,  [cx,go,~] = try_catch_load([new_dir_name 'centroid_x_' clus_name run_name],1,1);
      if go,  [cy,go,~] = try_catch_load([new_dir_name 'centroid_y_' clus_name run_name],1,1);
      if go, 
        num_clusters = nc.num_clusters; clear nc
        centroid_x = cx.centroid_x; clear cx
        centroid_y = cy.centroid_y; clear cy

        ngen = length(num_clusters);
        leng = ngen;
        lsx = (((basic_map_size(2)*2) -1)*2) -1;
        lsy = (((basic_map_size(1)*2) -1)*2) -1;

        u = 0;  v = 0;
        for gen = 1:ngen, %for each generation of this run
          u = v +1;
          v = sum(num_clusters(1:gen));
          refs = u:v;
          nc_dist = zeros(num_clusters(gen),1);
          x = centroid_x(refs); y = centroid_y(refs);
          for C = 1:num_clusters(gen),  
            distances = sqrt(((x(C)-x(:)).^2 +((y(C)-y(:)).^2)));
            [sd] = sort(distances);
            if num_clusters(gen)>1, 
              nc_dist(C) = sd(2);
            else, 
              nc_dist(C) = [];
              leng = leng -1;
            end
          end
          if length(nc_dist)>0, 
            m_nc_dist(gen,run) = mean(nc_dist);
          else, 
            m_nc_dist(gen,run) = 0;
          end
        end
        [tmncd] = find(m_nc_dist(1:ngen,r));
        this_m_nc_dist = m_nc_dist(tmncd,r);
        avg_nc_dist(b,r) = mean(this_m_nc_dist);

        area = lsx*lsy;
        rho = num_clusters(tmncd)'./area;
        mrd = 0.5*(rho.^-0.5);
        R = this_m_nc_dist./mrd;
        sigma_mrd = 0.26136./sqrt(num_clusters(tmncd)'.*rho);
        c = (this_m_nc_dist -mrd)./sigma_mrd;
        avgR(b,run) = mean(R);
        avgc(b,run) = mean(c);
        GEN_R(1:length(R),b) = R;
        GEN_c(1:length(c),b) = c;
      end
      end
      end
    end
  end
  only_these = find(avgR(b,:));
  AVG_R(b) = mean(avgR(b,only_these)); %simulations averaged of time averaged R
  STD_R(b) = std(avgR(b,only_these));  %simulations standard deviation of time averaged R
  AVG_c(b) = mean(avgc(b,only_these)); %simulations averaged of time averaged c
  STD_c(b) = std(avgc(b,only_these));  %simulations standard deviation of time averaged c
  %%
  if do_nearest_neighbors_plot==1
    if length(find(seed_distances(:,1)>max_size))>0
      par_dis = nearest_neighbor_distances; clear nearest_neighbor_distances
      len = lengths; clear lengths
      lengths = [len:step:(max(seed_distances(:,1))+step)];
      max_size = max(lengths);
      nearest_neighbor_distances = zeros(length(mutability),length(lengths));
      nearest_neighbor_distances(i,:) = hist(seed_distances(:,1),lengths);
      if i>1
        nearest_neighbor_distances(1:i-1,:) = [par_dis(1:i-1,:) zeros(length(1:i-1),(length(lengths)-length(len)))];
      end
    else
      nearest_neighbor_distances(i,:) = hist(seed_distances(:,1),lengths);
    end
    i = i +1;
    figure, semilogy(lengths,nearest_neighbor_distances(i-1,:),'x');
  end
end
end
end

AVG_NN_DIST = mean(avg_nn_dist,2);
STD_NN_DIST = std(avg_nn_dist')';
%%
if do_parent_distances_3d==1
  figure, mesh(log10(nearest_neighbor_distances));
  xlabel('100*(parent\_distances - overpop)');  ylabel('\mu');
  title(['parent\_distances\_' ...
    make_title_name(generalize_base_name(base_name),'')]);
end
%%
if do_dist_plot==1
%   load(['m_nn_dist_' base_name(1:length(base_name)-1)]);
  run_name = int2str(run);
  load(['population_' base_name run_name]);
  lsx = (((basic_map_size(1)*2)-1)*2)-1;
  lsy = (((basic_map_size(2)*2)-1)*2)-1;
  rho = population/(lsx*lsy);
  divide = 0.5*(rho.^-0.5);
  R = m_nn_dist((m_nn_dist(:,run)~=0),run)./divide';
  figure(3990157), plot(R,'x');
  avgR(run) = mean(R);
  u = 0;  v = 0;
  for gen = 1:ngen
    u = v +1; v = sum(population(1:gen));
    figure(3990157), loglog(population(gen),m_nn_dist(u:v,run),'x');  hold on;
  end
%     loglog(population,divide);
  xlabel('population');  ylabel('<parent distances>');
  title(make_title_name(base_name,''));
end
%%
if do_R_plot==1
  figure(3766);
  errorbar(mutability,AVG_R,STD_R,'x');
  title(make_title_name(make_data_name('',base_name,'',1),''));
  xlabel('\mu');  ylabel('R');
  
  figure(3767);
  errorbar(mutability,AVG_c,STD_c,'x');
  title(make_title_name(make_data_name('',base_name,'',1),''));
  xlabel('\mu');  ylabel('R significance for two-tailed test');
  hold on;
  plot([min(mutability) max(mutability)],[1.96 1.96],'k');
  plot([min(mutability) max(mutability)],[2.58 2.58],'b');
  plot([min(mutability) max(mutability)],[-1.96 -1.96],'k');
  plot([min(mutability) max(mutability)],[-2.58 -2.58],'b');
end