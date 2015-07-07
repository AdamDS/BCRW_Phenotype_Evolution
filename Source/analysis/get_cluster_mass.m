%% get_cluster_mass.m
% 
% -ADS 9*11*13
global SIMOPTS;
this_script = 'get_cluster_mass';
limit = SIMOPTS.limit;
tic;
N = numel(overpop)*numel(death_max)*numel(mutability);
AVG_S = zeros(N,1);
STD_S = zeros(N,1);
i = 0;
%%
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir,do_cd);
  [clus_name] = cluster_name(base_name);
  
  avg_s = zeros(length(SIMS),1);
  for run = SIMS, 
    run_name = int2str(run);

    fprintf('Attempting %s for %s \n',this_script,[clus_name run_name]);
    
    new_dir_name = split_cd(dir_name,run,split,make_dir,do_cd);

    r = find(run==SIMS);

    [p,go] = try_catch_load([new_dir_name 'population_' base_name run_name],1,1);
    if go,  [tx,go] = try_catch_load([new_dir_name 'trace_x_' base_name run_name],1,1);
    if go,  [ty,go] = try_catch_load([new_dir_name 'trace_y_' base_name run_name],1,1);
    if go,  [nc,go] = try_catch_load([new_dir_name 'num_clusters_' clus_name run_name],1,1);
    if go,  [tc,go] = try_catch_load([new_dir_name 'trace_cluster_' clus_name run_name],1,1);
    if go,  [onc,go] = try_catch_load([new_dir_name 'orgsnclusters_' clus_name run_name],1,1);
    if go,  
      population = p.population;  clear p
      trace_x = tx.trace_x; clear tx
      trace_y = ty.trace_y; clear ty
      num_clusters = nc.num_clusters; clear nc
      trace_cluster = tc.trace_cluster; clear tc
      orgsnclusters = onc.orgsnclusters;  clear onc

      finite = get_finite_clusters(population,num_clusters,trace_cluster,trace_x,trace_y,INFRAT);
      clear trace_x trace_y trace_cluster
      
      ngen = length(find(population>=limit));
      monc = max(orgsnclusters);
      range = limit:monc;
      s = zeros(ngen,1);

      for gen = 1:ngen, 
        [cu,cv] = gen_index(num_clusters,gen);
        f = finite(cu:cv); %f==1 is finite, f==0 is infinite along 1 axis, f==-1 is infinite both axes
        orgs = orgsnclusters(cu:cv);
        uo = unique(orgs(f==1));
        n = hist(orgs(f==1),range); %frequencies of finite s-clusters
        ns = n;%./population(gen); %cluster number distribution

        sns = range(n~=0).*ns(n~=0); %ns*s^2
        s(gen) = sum(range(n~=0).*sns)./sum(sns);
      end %for gen
      use = find(s~=0);
      avg_s(r) = mean(s(use));
    end %orgsnclusters
    end %trace_cluster
    end %num_clusters
    end %trace_y
    end %trace_x
    end %population
  end %for SIMS
  i = i +1;
  use = find(avg_s~=0);
  if ~isempty(use), 
    AVG_S(i) = mean(avg_s(use));
    STD_S(i) = std(avg_s(use));
  end
end %mu
end %dm
end %op
toc;

use = find(AVG_S~=0);
figure(3875); errorbar(mutability(use),AVG_S(use),STD_S(use));
xlabel('\mu');  ylabel('<cluster mass>');
xlim([min(mutability) max(mutability)]);
title(make_title_name(generalize_base_name(clus_name),''));

%% Old code
% global SIMOPTS;
% limit = SIMOPTS.limit;
% c = 'krbgkrbgkrbg';
% i = 0;
% N = numel(mutability)*numel(death_max)*numel(overpop);
% FISHER = zeros(1,N);
% SIGM = zeros(1,N);
% tic;
% for op = overpop, SIMOPTS.op = op;
% for dm = death_max, SIMOPTS.dm = dm;
% for mu = mutability, SIMOPTS.mu = mu;
%   make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
%   [clus_name] = cluster_name(base_name);
%   loc = zeros(numel(SIMS),1);
%   fre = [];
%   i = i +1;
%   r = 0;
%   go = ones(numel(SIMS),1);
%   [cnd,go] = try_catch_load(['cluster_number_distribution_' clus_name(1:(end-1))],1,1);
%   if go,  
%     cluster_number_distribution = cnd.cluster_number_distribution;  clear cnd,
%     S = [limit:(length(cluster_number_distribution)+limit-1)];
%     figure(6234); loglog(S,cluster_number_distribution,c(i)); hold on;
% %   for run = SIMS, 
% %     r = r +1;
% %     run_name = int2str(run);
% %     [onc,go(r)] = try_catch_load(['orgsnclusters_' clus_name run_name],go(r),1);
% %     if go(r), 
% %     orgsnclusters = onc.orgsnclusters;  clear onc
% %     [b,a] = hist(orgsnclusters,[limit:max(orgsnclusters)]);
% %     fre = cat_row(fre,b);
% %     loc(r) = max(orgsnclusters);
% %     end %orgsnclusters
% % %     end %population
% %   end %for sims
% %   if sum(go), 
% %   FRE = mean(fre);
% %   LOC = [limit:1:max(loc)];
% %   logFRE = log(FRE);
% %   logLOC = log(LOC);
% %   ok = find(FRE & LOC>ceil(0.1*max(loc)));%& FRE>1);%  %get asymptotic tail
% %   [fisher,intercept,sigm,sigb] = linear_fit(logLOC(ok),logFRE(ok),[],0);
% %   FISHER(i) = fisher;
% %   SIGM(i) = sigm;
% %   figure(6234); loglog(LOC,FRE,c(i));  hold on;
% %   end
%   end
% end %for mu
% end %for dm
% end %for op
% xlim([limit 100]);
% % if do_fisher_plot,  
% %   figure(77200 +10*reproduction +limit);
% %   if numel(mutability)>numel(death_max),  
% %     x_ = mutability;
% %   else, 
% %     x_ = death_max;
% %   end
% %   errorbar(x_,FISHER,SIGM);
% % end
% % cd 'G:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\Data'
% % load('G:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\Data\orgsnclusters_Bacterial_Mu_38_Uniform_2_Flatscape_1.mat')
% % oc1 = orgsnclusters;
% % load('G:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\Data\orgsnclusters_Bacterial_Mu_38_Uniform_2_Flatscape_2.mat')
% % oc2 = orgsnclusters;
% % load('G:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\Data\orgsnclusters_Bacterial_Mu_38_Uniform_2_Flatscape_3.mat')
% % oc3 = orgsnclusters;
% % load('G:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\Data\orgsnclusters_Bacterial_Mu_38_Uniform_2_Flatscape_4.mat')
% % oc4 = orgsnclusters;
% % load('G:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\Data\orgsnclusters_Bacterial_Mu_38_Uniform_2_Flatscape_5.mat')
% % oc5 = orgsnclusters;
% % [b1,a1] = hist(orgsnclusters(100));
% % [b1,a1] = hist(oc1(100));
% % [b2,a2] = hist(oc2(100));
% % [b3,a3] = hist(oc3(100));
% % [b4,a4] = hist(oc4(100));
% % [b5,a5] = hist(oc5(100));
% % [b1,a1] = hist(oc1(100),[3:2000]);
% % plot(at,b1,'x')