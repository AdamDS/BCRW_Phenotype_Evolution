%% get_nu_parallel.m *******************************************************
% Uses data collapse.
global SIMOPTS;
m = 0; c = 'ybrmgk'; c = [c c c c c];

num_sims = zeros(length(mutability),1);
sim_exists = lmlS;

win_pops = zeros(length(mutability),NGEN-1);

N = numel(mutability)*numel(death_max)*numel(overpop);

for bms = basic_map_sizes,  
  m = m +1;
  SIMOPTS.basic_map_size = [bms bms]; 
  maxpop = get_max_population();
  [A,xl,yl] = landscape_measures(SIMOPTS.basic_map_size);

  lmN = zeros(length(N),NGEN);
  AVG_GEN_POPS = lmN;
  STD_GEN_POPS = lmN;

  num_sims = zeros(N,1);
  sim_exists = lmlS;

  clear lmN

  alphas = zeros(N,1);
  sigas = zeros(N,1);
  
  i = 0;
  tic;
  for op = overpop, SIMOPTS.op = op;
  for dm = death_max, SIMOPTS.dm = dm;
  for mu = mutability, SIMOPTS.mu = mu;
    [base_name,dir_name] = NameAndCD(make_dir,do_cd);
    i = i +1;
    r = 0;
    %% from only_lt
    [p,go] = try_catch_load([dir_name 'populations_' base_name SIMrange],1,1);
    if go,  [lt,go] = try_catch_load([dir_name 'lifetimes_' base_name SIMrange],1,1);
    if go,  
      populations = p.populations;  clear p
      lifetimes = lt.lifetimes; clear lt
      gens = [ones(lenght(lifetimes),1) lifetimes];
      AVG_POPS(i) = mean(;
      AVG_GEN_POPS(i,:) = mean(populations,1);
      ngens(i,:) = lifetimes;
    end
    end
    r = length(lifetimes);
    %% from regular sims
    %build the simulation names for which you have chosen
    populations = zeros(length(SIMS),NGEN);
    for run = SIMS, 
      r = r +1;
      new_dir_name = split_cd(dir_name,run,split,make_dir,do_cd);
      exp_name = [base_name int2str(run)];
      pop_name = [new_dir_name 'population_' exp_name];
      [p,go] = try_catch_load(pop_name,1,1);
      if go, 
        num_sims(i) = num_sims(i) +1;
        sim_exists(i,r) = 1; %track which simulation data is used
        population = p.population;  clear p
        if ~loaded || no_bio, g = find(population>=limit);
        else, g = find(population); end
        populations(run,1:length(g)) = population(g);
        ngens(i,r) = length(g);
      end
      use = find(sim_exists(i,:)); %get sims that exist
      populations(~use,:) = []; %clear rows that have no data
      AVG_GEN_POPS(i,:) = mean(populations,1);
    end
    
    %% determine alpha
    RHO{i} = AVG_GEN_POPS(i,AVG_GEN_POPS(i,:)~=0)./maxpop;
    fracsur = length(find(ngens(i,:)==NGEN))./size(AVG_GEN_POPS,1); %fraction survived
    tmax = max(ngens(i,:));
    tmin = ceil(tmax/10);
    [alphas(i),~,sigas(i),~] = linear_fit(log10(tmin:tmax),log10(AVG_GEN_POPS(i,tmin:tmax)),[],0);
    fprintf('mu %1.4f\t alpha %1.4f\t +/- %1.4f\t fraction survived %1.2f\n',...
      mu,alphas(i),sigas(i),fracsur);
  end %mu
  end %dm
  end %op

  %% Determine critical point
  muc = input('What is the critical mutability (want alpha ~ 0.451)?'); %determine critical mu
  Delta = (mutability-muc);
  lessthan = find(Delta<0);
  greaterthan = find(Delta>0);
  alpha = alphas(min(abs(Delta)));
  
  %% Determine correlation time exponent: rho*t^alpha = f(t*Delta^nu_para)
  dt = 1;
  t = 1:dt:NGEN;  %1xNGEN
  dnp = 0.01;
  np = 1.10:dnp:1.50; %expect 1.295 ~ 1.30
  nnp = length(np);
  x = cell(N,nnp);
  y = cell(N,nnp);
  for j = 1:nnp,  %for each nu_para
    nu_para = np(j);
    for i = 1:N,  %for each time series
      if i~=find(Delta==0), %no need for collapsing on critical
        x{i}{j} = t(RHO{i}~=0)*Delta(i).^nu_para;
        y{i}{j} = RHO{i}.*t(RHO{i}~=0).^alpha;
      end
    end
  end
  
  % Bhatarcharjee & Seno 2001: goodness of data-collapse
  q = 1;  
  eta = 0.1;  %10% level of error
  % below critical
  figure(9011); hold on;
  I = length(lessthan); J = 1:(length(lessthan)-1);
  [npb(m),dnpb(m)] = goodness_of_data_collapse(np,x{lessthan}{:},y{lessthan}{:},I,J,q,eta);
  % above critical
  I = N-1;  J = (3+length(lessthan)):N; 
  [npa(m),dnpa(m)] = goodness_of_data_collapse(np,x{greaterthan}{:},y{greaterthan}{:},I,J,q,eta);
  % npb = nu_para below;  npa = nu_para above
  
  %% Determine correlation length exponent: rho*t^alpha = f(t/L^z)
  dz = 0.01;
  z = 0.37:dz:0.77;  %expect 0.567 ~ 0.57
  nz = length(z);
  x = cell(N,nz);
  y = cell(N,nz);
  L = min([xl yl]);
  for j = 1:nz,  %for each nu_perp
    Z = z(j);
    for i = 1:N,  %for each time series
      if i~=find(Delta==0), %no need for collapsing on critical
        x{i}{j} = t(RHO{i}~=0)./L.^Z;
        y{i}{j} = RHO{i}.*t(RHO{i}~=0).^alpha;
      end
    end
  end
  
  % Bhatarcharjee & Seno 2001: goodness of data-collapse
  q = 1;  
  eta = 0.1;  %10% level of error
  figure(907);
  I = N-1; J = 1:N;
  [nra(m),dnra(m)] = goodness_of_data_collapse(z,x,y,I,J,q,eta);  
  % nr = nu_perp
  
  
  

  figure(16161); plot(use,RHO{i}(use),c(i)); hold on;
  [slope(i),yint(i),eslope(i),eyint(i)] = linear_fit(log10(t(use2)),log10(RHO{i}(use2)),[],0);
  plot_nu_parallel(RHO{i}(use),t(use),mucrit,alpha,Delta,nu_parallel,c(i));
%   end
end %bms

toc;
% title(['Bacterial\_Mu\_x\_Uniform\_2\_Flatscape\_1000000\_NGEN, \mu_{c} = ' num2str(mucrit)],'FontSize',16);

set(gca,'yscale','log','xscale','log');
xlabel('t\Delta^{\nu}','FontSize',14); ylabel('\rho(t)t^{\alpha}','FontSize',14);
% xlim([]); 

%% Write to xls
% col = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'; 
% for i=1:11, 
%   xlswrite(['RHO_' int2str(use(1)) '_' int2str(diff(use(1:2))) '_' int2str(use(end)) ...
%     '_mc_335_delta_462_nupar_1_' int2str(mod(nu_parallel,1)*100) '_mu_33_001_34'],...
%     RHO{i}',...
%     [col(i) int2str(2) ':' col(i) int2str(length(RHO{i}))]); 
% end


%% Old code: 
% 
% %% Alternative scaling
% % Tmax = 10^4;
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_33_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 1; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_331_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 2; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_332_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 3; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_333_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 4; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_334_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 5; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_335_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 6; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_336_Uniform_2_Flatscape_1000000_NGEN_1_100.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 7; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_337_Uniform_2_Flatscape_1000000_NGEN_1_100_300_359.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 8; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_338_Uniform_2_Flatscape_1000000_NGEN_300_359.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 9; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_339_Uniform_2_Flatscape_1000000_NGEN_300_359.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 10; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_34_Uniform_2_Flatscape_1000000_NGEN_300_365.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 11; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_341_Uniform_2_Flatscape_1000000_NGEN_305_319_325_329_355_359.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 12; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_342_Uniform_2_Flatscape_1000000_NGEN_305_319.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 13; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % load('H:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\1000000_NGEN\Data\lifetimes_Bacterial_Mu_343_Uniform_2_Flatscape_1000000_NGEN_310_314.mat');
% % for i=1:length(lifetimes), if lifetimes(i)>Tmax, lifetimes(i) = Tmax; end,  end
% % i = 14; 
% % T(i) = mean(lifetimes); varT(i) = std(lifetimes);
% % 
% % mu = 0.330:0.001:0.343;
% % mucrit = 0.335; eps = abs(mu(mu<mucrit)-mucrit);
% % figure(1), plot(eps,T(1:length(eps)),'-x');  set(gca,'yscale','log','xscale','log');
% % figure(2), plot(eps,varT(1:length(eps)),'-x'); set(gca,'yscale','log','xscale','log');
% 
% 
% %%
% if only_lt, 
%     
%   else, 
%     %build the simulation names for which you have chosen
%     gen_pops = zeros(length(SIMS),NGEN);
%     i = i +1;
%     r = 0;
%     for run = SIMS, 
%       r = r +1;
%       new_dir_name = split_cd(dir_name,run,split,make_dir,do_cd);
%       run_name = int2str(run);
%       exp_name = [base_name run_name];
%       pop_name = [new_dir_name 'population_' exp_name];
%       [p,go] = try_catch_load(pop_name,1,1);
%       if go, 
%         num_sims(i) = num_sims(i) +1;
%         sim_exists(i,r) = 1; %track which simulation data is used
%         population = p.population;  clear p
%         g = find(population>=limit);
%         actual_ngen = max(g);
%         ngens(i,r) = actual_ngen;
%         avg_pop(i,run) = mean(population(g));
%         gen_pops(run,1:length(g)) = population(g);
%       end
%     end
%   end
%
%
%   use = find(sim_exists(i,:)); these = 1:length(use);
%   gen_pops(~use,:) = [];
%     
%   [rgood(i),cgood] = size(gen_pops);
%   
%   AVG_POPS(i) = mean(avg_pop(i,:));
%   STD_POPS(i) = std(avg_pop(i,:));
% 
%   for gen=1:max(ngens(i,:)),  
%     alive = find(gen_pops(:,gen));
%     AVG_GEN_POPS(i,gen) = mean(gen_pops(alive,gen),1);
%     STD_GEN_POPS(i,gen) = std(gen_pops(alive,gen));
%   end
%   
% %   for gen = 1:(NGEN-window),  
% %     win_pops(run,gen) = mean(AVG_GEN_POPS(i,gen:(gen+window)),2);
% %   end
% %   for gen = 1:NGEN, 
% %     use = find(gen_pops(:,gen));
% %     if isempty(use),  gen = NGEN-1; 
% %     else, sur_pops(i,gen) = mean(gen_pops(use,gen));  end
% %   end
% 
% 
% %% LOAD STUFF
% % %   color = c(i); 
% % close all;
% % clear all;
% % clc
% % % C = 'kkkkkk';
% % C = 'brykg'; 
% % c = [C C C];
% % mucrit = 0.335;
% % delta = 0.462; %0.462 +- 0.033 from mu=0.335 for 100sims
% % nu_parallel = 1.19; %0.335 1.05
% % A = 37544;%45*45;%/37544;
% % t = 1:10^6;
% % use = 200:200:100000;
% % mutability = 0.33:0.001:0.34;
% % RHO = cell(numel(mutability),1);
% % T = t(use);
% % 
% % i = 1;
% % Delta = abs(mutability(i)-mucrit);
% % load('populations_Bacterial_Mu_33_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % rho = mean(populations)./A;
% % % use = find(rho);
% % RHO{i} = rho(use);
% % [slope(i),yint(i),eslope(i),eyint(i)] = linear_fit(log10(T),log10(RHO{i}),[],0);
% % plot_nu_parallel(RHO{i},T,mucrit,delta,Delta,nu_parallel,c(i));
% % 
% % i = 2;
% % Delta = abs(mutability(i)-mucrit);
% % load('populations_Bacterial_Mu_331_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % rho = mean(populations)./A;
% % % use = find(rho);
% % RHO{i} = rho(use);
% % [slope(i),yint(i),eslope(i),eyint(i)] = linear_fit(log10(T),log10(RHO{i}),[],0);
% % plot_nu_parallel(RHO{i},T,mucrit,delta,Delta,nu_parallel,c(i));
% % 
% % i = 3;
% % Delta = abs(mutability(i)-mucrit);
% % load('populations_Bacterial_Mu_332_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % rho = mean(populations)./A;
% % % use = find(rho);
% % RHO{i} = rho(use);
% % [slope(i),yint(i),eslope(i),eyint(i)] = linear_fit(log10(T),log10(RHO{i}),[],0);
% % plot_nu_parallel(RHO{i},T,mucrit,delta,Delta,nu_parallel,c(i));
% % 
% % i = 4;
% % Delta = abs(mutability(i)-mucrit);
% % load('populations_Bacterial_Mu_333_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % rho = mean(populations)./A;
% % % use = find(rho);
% % RHO{i} = rho(use);
% % [slope(i),yint(i),eslope(i),eyint(i)] = linear_fit(log10(T),log10(RHO{i}),[],0);
% % plot_nu_parallel(RHO{i},T,mucrit,delta,Delta,nu_parallel,c(i));
% % 
% % i = 5;
% % Delta = abs(mutability(i)-mucrit);
% % load('populations_Bacterial_Mu_334_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % rho = mean(populations)./A;
% % % use = find(rho);
% % RHO{i} = rho(use);
% % [slope(i),yint(i),eslope(i),eyint(i)] = linear_fit(log10(T),log10(RHO{i}),[],0);
% % plot_nu_parallel(RHO{i},T,mucrit,delta,Delta,nu_parallel,c(i));
% % 
% % i = 6;
% % Delta = abs(mutability(i)-mucrit);
% % load('populations_Bacterial_Mu_335_Uniform_2_Flatscape_1000000_NGEN_101_200.mat');
% % rho = mean(populations)./A;
% % % use = find(rho);
% % RHO{i} = rho(use);
% % [slope(i),yint(i),eslope(i),eyint(i)] = linear_fit(log10(T),log10(RHO{i}),[],0);
% % plot_nu_parallel(RHO{i},T,mucrit,delta,Delta,nu_parallel,c(i));
% % 
% % i = 7;
% % Delta = abs(mutability(i)-mucrit);
% % load('populations_Bacterial_Mu_336_Uniform_2_Flatscape_1000000_NGEN_1_100.mat');
% % rho = mean(populations)./A;
% % % use = find(rho);
% % RHO{i} = rho(use);
% % [slope(i),yint(i),eslope(i),eyint(i)] = linear_fit(log10(T),log10(RHO{i}),[],0);
% % plot_nu_parallel(RHO{i},T,mucrit,delta,Delta,nu_parallel,c(i));
% % 
% % i = 8;
% % Delta = abs(mutability(i)-mucrit);
% % load('populations_Bacterial_Mu_337_Uniform_2_Flatscape_1000000_NGEN_1_100_300_359.mat');
% % rho = mean(populations)./A;
% % % use = find(rho);
% % RHO{i} = rho(use);
% % [slope(i),yint(i),eslope(i),eyint(i)] = linear_fit(log10(T),log10(RHO{i}),[],0);
% % plot_nu_parallel(RHO{i},T,mucrit,delta,Delta,nu_parallel,c(i));
% % 
% % i = 9;
% % Delta = abs(mutability(i)-mucrit);
% % load('populations_Bacterial_Mu_338_Uniform_2_Flatscape_1000000_NGEN_300_359.mat');
% % rho = mean(populations)./A;
% % % use = find(rho);
% % RHO{i} = rho(use);
% % [slope(i),yint(i),eslope(i),eyint(i)] = linear_fit(log10(T),log10(RHO{i}),[],0);
% % plot_nu_parallel(RHO{i},T,mucrit,delta,Delta,nu_parallel,c(i));
% % 
% % i = 10;
% % Delta = abs(mutability(i)-mucrit);
% % load('populations_Bacterial_Mu_339_Uniform_2_Flatscape_1000000_NGEN_300_359.mat');
% % rho = mean(populations)./A;
% % % use = find(rho);
% % RHO{i} = rho(use);
% % [slope(i),yint(i),eslope(i),eyint(i)] = linear_fit(log10(T),log10(RHO{i}),[],0);
% % plot_nu_parallel(RHO{i},T,mucrit,delta,Delta,nu_parallel,c(i));
% % 
% % i = 11;
% % Delta = abs(mutability(i)-mucrit);
% % load('populations_Bacterial_Mu_34_Uniform_2_Flatscape_1000000_NGEN_300_365.mat');
% % rho = mean(populations)./A;
% % % use = find(rho);
% % RHO{i} = rho(use);
% % [slope(i),yint(i),eslope(i),eyint(i)] = linear_fit(log10(T),log10(RHO{i}),[],0);
% % plot_nu_parallel(RHO{i},T,mucrit,delta,Delta,nu_parallel,c(i));
% % 
% %%