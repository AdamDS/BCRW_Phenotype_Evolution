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

global SIMOPTS;
% max_size = 3;%(((max(basic_map_size)*2)-1)*2)-1;
AVG_GEN_NN_DIST = zeros(NGEN,length(mutability));
STD_GEN_NN_DIST = zeros(NGEN,length(mutability));
avg_nn_dist = zeros(length(mutability),length(SIMS));
AVG_NN_DIST = zeros(length(mutability),1);
STD_NN_DIST = zeros(length(mutability),1);
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
    exp_name = [base_name run_name];
    pn = ['population_' exp_name];  load(pn);
    skip = 1;
    try, 
      load(['seed_distances_' exp_name]), 
      skip = 0;
    catch esd, 
      esd.message
      skip = 1;
    end
    if skip==0, 
      tcsn = ['trace_cluster_seed_' exp_name];  load(tcsn);
      txn = ['trace_x_' exp_name];  load(txn);
      tyn = ['trace_y_' exp_name];  load(tyn);
      
      ngen = length(population);

      hash_length = mu;
      lsx = (((basic_map_size(2)*2) -1)*2) -1;
      lsy = (((basic_map_size(1)*2) -1)*2) -1;
      hlmx = mod(lsx,hash_length);
      hlmy = mod(lsy,hash_length);
      grid_horizontal = [0:hash_length:lsx+hlmx];
      grid_vertical = [0:hash_length:lsy+hlmy];
      areas = zeros(1,ngen);

      u = 0;  v = 0;
      for gen = 1:ngen
        u = v +1;
        v = sum(population(1:gen));
        if population(gen)>=limit
          refs = u:v;
          if long==1
            mates = u +double(trace_cluster_seed(refs,1)) -1;
            seed_distances(refs,1) = sqrt((trace_x(refs) -trace_x(mates)).^2 ...
              +(trace_y(refs) -trace_y(mates)).^2);
            alternates = u +double(trace_cluster_seed(refs,2)) -1;
            seed_distances(refs,2) = sqrt((trace_x(refs)-trace_x(alternates)).^2 ...
              +(trace_y(refs) -trace_y(alternates)).^2);
          end
          m_nn_dist(gen,run) = mean(seed_distances(refs,1));
        end
      end
      if long==1, save(['seed_distances_' exp_name],'seed_distances');  end
      mnd = m_nn_dist(find(m_nn_dist(:,run)),run);

      this_nn_dist = m_nn_dist(find(m_nn_dist(:,r)),r);
      avg_nn_dist(b,r) = mean(this_nn_dist);
      %%
      areas = lsx*lsy*ones(size(areas));
      rho = population./areas;
      mrd = 0.5*(rho.^-0.5)';
      R = this_nn_dist./mrd;
      sigma_mrd = 0.26136./sqrt(population.*rho)';
      c = (this_nn_dist -mrd)./sigma_mrd;
      avgR(b,run) = mean(R);
      avgc(b,run) = mean(c);
      if long==0
        GEN_R(1:length(R),b) = R;
        GEN_c(1:length(c),b) = c;
      else
        GEN_R(1:length(R),b) = R';
        GEN_c(1:length(c),b) = c';
      end
      AVG_GEN_NN_DIST(:,b) = mean(m_nn_dist,2);
      STD_GEN_NN_DIST(:,b) = std(m_nn_dist')';
      %%
      AVG_R(b) = mean(avgR(b,:));
      STD_R(b) = std(avgR(b,:));
      AVG_c(b) = mean(avgc(b,:));
      STD_c(b) = std(avgc(b,:));
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
%   save(['m_nn_dist_' exp_name],'m_nn_dist');
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