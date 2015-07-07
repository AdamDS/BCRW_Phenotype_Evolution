global SIMOPTS;
%% get_kills
avg_kills = zeros(length(mutability),3,length(SIMS));
AVG_KILLS = zeros(length(mutability),3);
STD_KILLS = zeros(length(mutability),3);
avg_riv = zeros(length(mutability),length(SIMS));
AVG_RIV = zeros(length(mutability),3);
STD_RIV = zeros(length(mutability),3);
std_op = zeros(length(mutability),length(SIMS));
std_rd = zeros(length(mutability),length(SIMS));
std_cj = zeros(length(mutability),length(SIMS));

avg_pk = zeros(length(mutability),3,length(SIMS));
AVG_PK = zeros(length(mutability),3);
STD_PK = zeros(length(mutability),3);
avg_pr = zeros(length(mutability),length(SIMS));
AVG_PR = zeros(length(mutability),3);
STD_PR = zeros(length(mutability),3);
std_pop = zeros(length(mutability),length(SIMS));
std_prd = zeros(length(mutability),length(SIMS));
std_pcj = zeros(length(mutability),length(SIMS));

b = 0;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  b = b +1;
  R = 0;
  for run = SIMS, 
    R = R +1;
    run_name = int2str(run);
    new_dir_name = split_cd(dir_name,run,split,make_dir,do_cd);
    exp_name = [base_name run_name];
    go = 1; [p,go] = try_catch_load([new_dir_name 'population_' exp_name],go,1);
  if go==1, [k,go] = try_catch_load([new_dir_name 'kills_' exp_name],go,1);
  if go==1, [r,go] = try_catch_load([new_dir_name 'rivalries_' exp_name],go,1);
  if go==1, 
    population = p.population;  clear p      
    kills = k.kills;  clear k
    ngen = length(population);
    if ngen<transience, start = floor(ngen/2);
    else, start = transience+1; end
    avg_kills(b,:,R) = mean(kills(start:ngen,:),1);
    rivalries = r.rivalries;  clear r
    avg_riv(b,R) = mean(rivalries(start:ngen),1);
    B = 2*population(start:ngen)';
    C = B-kills(start:ngen,2);
    A = C-kills(start:ngen,3);
    avg_pk(b,:,R) = mean([kills(start:ngen,1)./B ...
                          kills(start:ngen,2)./C ...
                          kills(start:ngen,3)./A],1);
    avg_pr(b,R) = mean(rivalries(start:ngen)./B,1);
  end
  end
  end
  end
end
end
end
AVG_KILLS = mean(avg_kills,3);
std_op(:,SIMS) = avg_kills(:,1,:);  std_op = std(std_op')';
std_rd(:,SIMS) = avg_kills(:,2,:);  std_rd = std(std_rd')';
std_cj(:,SIMS) = avg_kills(:,3,:);  std_cj = std(std_cj')';
STD_KILLS = [std_op std_rd std_cj]; clear std_op std_rd std_cj;
AVG_RIV = mean(avg_riv,2);
STD_RIV = std(avg_riv')';

AVG_PK = mean(avg_pk,3);
std_pop(:,SIMS) = avg_pk(:,1,:);  std_pop = std(std_pop')';
std_prd(:,SIMS) = avg_pk(:,2,:);  std_prd = std(std_prd')';
std_pcj(:,SIMS) = avg_pk(:,3,:);  std_pcj = std(std_pcj')';
STD_PK = [std_pop std_prd std_pcj]; clear std_pop std_prd std_pcj;
AVG_PR = mean(avg_pr,2);
STD_PR = std(avg_pr')';
if do_kills_plot==1
  figure(666); 
  errorbar(mutability,AVG_KILLS(:,1),STD_KILLS(:,1),'x'); hold on;
  if all_kills==1, 
    errorbar(mutability,AVG_KILLS(:,2),STD_KILLS(:,2),'d');
    errorbar(mutability,AVG_KILLS(:,3),STD_KILLS(:,3),'s');
  end
  errorbar(mutability,AVG_RIV,STD_RIV,'*');
  xlabel('\mu');  ylabel('kills');  
  if all_kills==1, 
    legend('Total Competition','Random','Wanderers','Sibling Rivalry');
  else, 
    legend('Total Competition','Sibling Rivalry');
  end
  title(['average kills for ' ...
    make_title_name(generalize_base_name(base_name),'')]);
end
if do_kills_plot==1
  figure(999); 
  errorbar(mutability,AVG_PK(:,1),STD_PK(:,1),'x'); hold on;
  if all_kills==1, 
    errorbar(mutability,AVG_PK(:,2),STD_PK(:,2),'d');
    errorbar(mutability,AVG_PK(:,3),STD_PK(:,3),'s');
  end
  errorbar(mutability,AVG_PR,STD_PR,'*');
  xlabel('\mu');  ylabel('kills');  
  if all_kills==1, 
    legend('Total Competition','Random','Wanderers','Sibling Rivalry');
  else, 
    legend('Total Competition','Sibling Rivalry');
  end
  title(['average kills for ' ...
    make_title_name(generalize_base_name(base_name),'')]);
end