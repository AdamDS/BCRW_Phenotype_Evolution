%% get_kills
global SIMOPTS;
avg_kills = zeros(length(mutability),3,length(SIMS));
AVG_KILLS = zeros(length(mutability),3);
STD_KILLS = zeros(length(mutability),3);
avg_riv = zeros(length(mutability),length(SIMS));
AVG_RIV = zeros(length(mutability),3);
STD_RIV = zeros(length(mutability),3);
% SR_NR_POP = zeros(length(mutability),3,length(SIMS));

b = 0;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  b = b +1;
  %build the simulation names for which you have chosen to generate cluster data
  r = 0;
  for run = SIMS
    r = r +1;
    run_name = int2str(run);
    exp_name = [base_name run_name];
    pn = ['population_' exp_name];  
    if exist([pn '.mat'])==2
      load(pn);
%       fix_population(base_name,run_name,limit);
      load(pn);
      kn = ['kills_' exp_name];  load(kn);

      rn = ['rivalries_' exp_name];  load(rn);
      avg_riv(b,r) = mean(rivalries,1);
      neighbor_rivalries = kills(:,1) -rivalries;
      sr_nr_pop = [rivalries neighbor_rivalries population'];
%       SR_NR_POP(:,:,run) = sr_nr_pop;      
      figure(5181199);
      plot(population,rivalries,'x'); hold on;
      figure(9314802);
      plot(population,neighbor_rivalries,'x'); hold on;
      pause;
    else
      ['Does not exist: ' base_name run_name]
    end
  end
end
end
end
% AVG_KILLS = mean(avg_kills,3);
% std_op(:,SIMS) = avg_kills(:,1,:);  std_op = std(std_op')';
% std_rd(:,SIMS) = avg_kills(:,2,:);  std_rd = std(std_rd')';
% std_cj(:,SIMS) = avg_kills(:,3,:);  std_cj = std(std_cj')';
% STD_KILLS = [std_op std_rd std_cj]; clear std_op std_rd std_cj;
% AVG_RIV = mean(avg_riv,2);
% STD_RIV = std(avg_riv')';
% if do_kills_plot==1
%   figure(666); 
%   errorbar(mutability,AVG_KILLS(:,1),STD_KILLS(:,1),'x'); hold on;
% %   errorbar(mutability,AVG_KILLS(:,2),STD_KILLS(:,2),'d');
% %   errorbar(mutability,AVG_KILLS(:,3),STD_KILLS(:,3),'s');
%   errorbar(mutability,AVG_RIV,STD_RIV,'*');
%   xlabel('\mu');  ylabel('kills');  %legend('Total Competition','Random','Wanderers','Sibling Rivalry');
%   legend('Total Competition','Sibling Rivalry');
%   title(['average kills for ' ...
%     make_title_name(generalize_base_name(base_name),'')]);
% end