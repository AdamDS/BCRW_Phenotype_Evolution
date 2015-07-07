%% get_population_decay_data.m (script)
lmN = zeros(length(N),NGEN);
AVG_GEN_POPS = zeros(1,NGEN);
STD_GEN_POPS = lmN;

num_sims = zeros(N,1);
sim_exists = zeros(N,numel(SIMS));
ngens = zeros(N,numel(SIMS));

clear lmN

T = zeros(N,2);

RHO = cell(N,1);

fprintf('\t\tmu\t\t\t\t alpha\t\t\t\t +/-\t\t survived\t num_sims\n');
i = 0;
tic;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  if do_default_density,  %IPOP above assumed as if for 12x12
    IPOP = scale_IPOP(); %IPOP = IPOP/37544*(max_population(bms))
  end

  [base_name,dir_name] = NameAndCD(make_dir,do_cd);

  i = i +1;
  r = 0;
  %% from only_lt
  lifetimes = [];
  [p,gog] = try_catch_load([dir_name 'populations_' base_name SIMrange],1,0);
  if gog,  [lt,gog] = try_catch_load([dir_name 'lifetimes_' base_name SIMrange],1,0);
  if gog,  
    populations = p.populations;  clear p
    lifetimes = lt.lifetimes; clear lt
    gens = [ones(length(lifetimes),1) lifetimes];
    AVG_GEN_POPS = mean(populations,1);
    num_sims(i) = num_sims(i) +length(lifetimes);
    r = num_sims(i);
    ngens(i,1:num_sims(i)) = lifetimes;
    clear populations lifetimes
  end
  end

  %% from regular sims
  if num_sims(i)~=length(SIMS), 
    %build the simulation names for which you have chosen
    populations = zeros(length(SIMS),NGEN);
    for run = SIMS, 
      r = r +1;
      new_dir_name = split_cd(dir_name,run,split,make_dir,do_cd);
      exp_name = [base_name int2str(run)];
      pop_name = [new_dir_name 'population_' exp_name];
      [p,goi] = try_catch_load(pop_name,1,0);
      if goi, 
        num_sims(i) = num_sims(i) +1;
        sim_exists(i,r) = 1; %track which simulation data is used
        population = p.population;  clear p
        if ~loaded || no_bio, g = find(population>=limit);
        else, g = find(population); end
        populations(run,1:length(g)) = population(g);
        ngens(i,r) = length(g);
      end
      use = find(sim_exists(i,:)); %get sims that exist
      if length(use)>0, 
        populations(~use,:) = []; %clear rows that have no data
        AVG_GEN_POPS = mean(populations,1);
      end
    end
  end
  mlt(i,m) = mean(ngens(i,ngens(i,:)~=0));
  %% determine alpha
  RHO{i} = AVG_GEN_POPS(AVG_GEN_POPS(1:NGEN)~=0)./maxpop;
  fracsur(i,m) = length(find(ngens(i,:)>=NGEN))./num_sims(i); %fraction survived
  fT(i,:) = [amM(1) amM(2)];  %alpha linear fit times  
%   T(i,2) = fT(i,2);%max(ngens(i,:));%fT(i,2);%
%   T(i,1) = ceil(T(i,2)/10);%fT(i,1);%ceil(T(i,2)/10);%fT(i,1);%
  denfig = 1620+m;
  figure(denfig); hold on;
%   
  samp = 1:length(RHO{i});
  plot(log10(samp),log10(RHO{i}(samp)),c(i));
%   plot(samp,RHO{i},c(i));
%   
  if fT(i,2)>length(RHO{i}), fT(i,2) = length(RHO{i}); end
  samp = fT(i,1):dt:fT(i,2);
  [alpha_decay(i,m),a_int(i,m),sigas(i,m),~,chi2(i,m)] = linear_fit( ...
    log10(samp),log10(RHO{i}(samp)),[],0);%denfig);

  fprintf('\t\t%1.3f\t\t\t %1.4f\t\t\t %1.4f\t\t\t %1.2f\t\t\t %1.0f\n', ...
    mutability(i),-alpha_decay(i,m),sigas(i,m),fracsur(i,m),num_sims(i));
end %mu
end %dm
end %op
fprintf('\n');

alphaxls = ['decay_' generalize_base_name(base_name)];
T = [1:10000]';
lT = length(T)+1;
xlswrite(alphaxls,T,'rawN',['A2:A' int2str(lT)]);
alphabet = 'BCDEFGHIJKLMN';
for i=1:N,  
  T = [1:10000]';
  lRHO = length(RHO{i});
  if lRHO<T(end), T(length(RHO{i})+1:end) = []; finxls = 1+lRHO;
  else, finxls = 1+T(end);  end
  xlswrite(alphaxls,RHO{i}(T)','rawN',[alphabet(i) '2:' alphabet(i) int2str(finxls)]);
  xlswrite(alphaxls,[mutability' a_int(:,m) alpha_decay(:,m)],'fit_params');
end

[~,crit] = min(abs(mutability-muc));
plot(log10(samp),(a_int(crit,m)+alpha_decay(crit,m)*log10(samp)),'k');

set(gca,'xscale','log','yscale','log');
title([int2str(xl(m)) 'x' int2str(yl(m)) ' densities']);
xlabel('Generations');  ylabel('\rho');