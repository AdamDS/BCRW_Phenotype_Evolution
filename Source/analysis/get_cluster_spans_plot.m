%% get_cluster_spans_plot.m
% 
% -ADS 9*25*13
global SIMOPTS;
this_script = 'get_cluster_spans_plot';
limit = SIMOPTS.limit;
tic;
N = numel(overpop)*numel(death_max)*numel(mutability);
i = 0;
%%
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  [base_name,dir_name] = NameAndCD(make_dir,do_cd);
  [clus_name] = cluster_name(base_name);
  
  for run = SIMS, 
    run_name = int2str(run);
    new_dir_name = split_cd(dir_name,run,split,make_dir,do_cd);
    r = find(run==SIMS);
    fprintf('Attempting %s for %s \n',this_script,[clus_name run_name]);
    
    [cs,go] = try_catch_load([new_dir_name 'cluster_spans_' clus_name run_name],1,1);
    if go,  [onc,go] = try_catch_load([new_dir_name 'orgsnclusters_' clus_name run_name],1,1);
    if go,  [gr,go] = try_catch_load([new_dir_name 'gyration_radii_' clus_name run_name],1,1);
    if go,  
      cluster_spans = cs.cluster_spans; clear cs
      orgsnclusters = onc.orgsnclusters;  clear onc
      gyration_radii = gr.gyration_radii; clear gr

      figure(10000*mu+100000*op+100*dm);
%       subplot(3,1,1); 
        plot(orgsnclusters,cluster_spans(:,1),'xk'); hold on;
        plot(orgsnclusters,cluster_spans(:,2),'xb'); hold off;
%       subplot(3,1,2);
%         plot(orgsnclusters,gyration_radii);
%       subplot(3,1,3);
%         plot(gyration_radii,cluster_spans(:,1),'k');  hold on;
%         plot(gyration_radii,cluster_spans(:,2),'b');  hold off;
    end %orgsnclusters
    end %trace_cluster
    end %num_clusters
    
  end %for SIMS
%   i = i +1;
%   use = find(avg_s~=0);
%   if ~isempty(use), 
%     AVG_S(i) = mean(avg_s(use));
%     STD_S(i) = std(avg_s(use));
%   end
end %mu
end %dm
end %op
toc;

% use = find(AVG_S~=0);
% figure(3875); errorbar(mutability(use),AVG_S(use),STD_S(use));
% xlabel('\mu');  ylabel('<cluster mass>');
% xlim([min(mutability) max(mutability)]);
% title(make_title_name(generalize_base_name(clus_name),''));

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