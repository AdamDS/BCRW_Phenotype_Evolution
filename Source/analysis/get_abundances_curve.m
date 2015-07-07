% function [abundances,num_bins] = get_abundances_curve(do_abundances_plot,log_abundance)
global SIMOPTS;
i = 0;  flag = 1;
AVG_AB = [];  avg_abs = [];
num_abundances = zeros(length(mutability),1);
ABUNDANCES = [];
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  % single simulation analysis
  abundance = []; abundances = [];
  for run = SIMS
    run_name = int2str(run);

    load(['num_clusters_' base_name run_name]);
    load(['orgsnclusters_' base_name run_name]);
    ngen = length(find(num_clusters));
%     these_onc = orgsnclusters(after_transience(1):length(orgsnclusters));
    num_bins = max(num_clusters);
    abundance = hist(orgsnclusters,num_bins);
    if do_abundances_single==1
      if do_abundances_plot==1, 
        figure(9870000+mu*1000+limit); 
        plot(1:num_bins,abundances,'x'); 
        title(make_title_name(base_name,run_name));
        xlabel('orgsnclusters'); ylabel('abundance');
      end
      if do_abundances_log==1
        figure(8760000+mu*1000+limit);
        semilogy(abundances,'x');
        title(make_title_name(base_name,run_name));
        xlabel('orgsnclusters'); ylabel('log_{10} abundance');
      end
    end
    la = length(abundance);
    las = size((abundances),2);
    if la>las
      if flag==0
        qwer = abundances;  clear abundances
        abundances = [[qwer zeros(size(qwer,1),la-las)]; abundance];
      else
        abundances = [[abundances]; abundance];
        flag = 0;
      end
    else
      abundances = [abundances; [abundance zeros(1,las-la)]]; %#ok<AGROW>
    end
  end
  i = i +1;
  avg_ab_bn = mean(abundances,1); %average over all runs
  std_ab_bn = std(abundances,1);  %std over all runs
  if do_abundances_3d==1
    laab = length(avg_ab_bn); %length of avg_ab_bn
    laa = size(avg_abs,2);    %length of combined avg_ab_bn in avg_abs
    if laa >= laab  %if the new averages vector is shorter, then append zeroes to it
      avg_abs = [avg_abs; avg_ab_bn zeros(1,laa-laab)];
      std_abs = [std_abs; std_ab_bn zeros(1,laa-laab)];
    else  %if the new averages vector is longer, then append zeroes to the combined averages vector
      avg_abs = [avg_ab_bn zeros(size(avg_ab_bn,1),size(avg_ab_bn,2)+(laab-laa))]; clear avg_ab_bn
      std_abs = [std_ab_bn zeros(size(std_ab_bn,1),size(std_ab_bn,2)+(laab-laa))]; clear std_ab_bn
    end
  end
%   nz = find(avg_ab_bn);
  AB = [3:length(avg_ab_bn)];
  AVG_AB = avg_ab_bn;
  num_abundances(i) = length(avg_ab_bn);
  ABUNDANCES = [ABUNDANCES avg_ab_bn]; %#ok<AGROW>
%   AVG_AB = avg_ab_bn(nz);
%   STD_AB = std_ab_bn(nz);
  % AVG_AB = avg_ab_bn;
  if do_abundances_3d==1,
    figure(990000+mu*1000+limit);
    surf(avg_abs);
  end
  if do_abundances_general==1
    nz = find(avg_ab_bn);
    AVG_AB = avg_ab_bn(nz);
    if do_abundances_plot==1
      figure(150000+mu*1000+limit);
      plot(nz,AVG_AB,'x');
      title(make_title_name(base_name,''));  
      xlabel('orgsnclusters');  ylabel('<abundance>');
    end
    if do_abundances_log_error==1
      figure(160000+mu*1000+limit);
      lAVG_AB = log10(AVG_AB);
      lSTD_AB = log10(STD_AB);
      errorbar(nz,lAVG_AB,lSTD_AB,'x');
      title(make_title_name(base_name,''));  
      xlabel('orgsnclusters');  ylabel('log_{10} <abundance>');
      xlim([-1 max(nz)+1]);
    end
    if do_abundances_log==1
      figure(170000+mu*1000+limit);
      semilogy(nz,AVG_AB,'x');
      title(make_title_name(base_name,''));  
      xlabel('orgsnclusters');  ylabel('<abundance>');
      xlim([-1 max(nz)+1]);
    end
    if do_abundances_avg==1
      figure(180000+mu*1000+limit);
      errorbar(nz,AVG_AB,std(abundances(:,nz)),'x');
      title(make_title_name(base_name,''));  
      xlabel('orgsnclusters');  ylabel('abundance');
    end
  end
end
end
end
if record_abundances==1
  save_abundances;
end